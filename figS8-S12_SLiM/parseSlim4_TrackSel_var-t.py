"""
Usage:
python %prog <slim_output> <parsed_filename> <scaling_lambda>

=======================================================
This script is specifically tailored to parse output from `Slim4_200kb_HC_neutral_ConstN1e4_t9x500gen_var-n_Lambda1.eidos` or `Slim4_200kb_HC_selection_ConstN1e4_t9x500gen_var-n_Lambda1.eidos`, so it assumes particular patterns in the output, and takes certain parameter values as given:
    - In the slim scripts for de novo sweeps, selection started at generation 199999 and lasted till 204000
    - Both scripts use lambda=1, so the population size is 1e4 diploid individuals
    - This script make slim records (output under "T") the number of copies of the selected allele at every generation since selection starts
    - Simulation restarts with another seed if this allele is lost (script writes out "Mutation lost")
    - Output will include the population allele count trajectory as a comment before the header
"""
import sys, re, time


# go thru the slim output
def readInput(infile):
    Pop = []  # list of dictionaries
    RefPop = {}
    Sizes = []  # sample sizes
    times = []  # sample times
    numSites = []  # for sanity check
    all_positions = []
    dup = []
    # for tracking allele freq
    last_id = '0'
    t0 = 0
    gen = []
    # counts = [0 for meh in range(4001)]  # list of copies per generation after introduction (from 199999 to 204000)
    counts = []

    restart_counter = 0
    if re.search(r'[N|n]eut', infile):
        neutral_flag = True
    else:
        neutral_flag = False

    target_pos = None
    line_number = 0

    # open file
    slimout = open(infile, 'r')
    l = slimout.readline()
    line_number += 1

    # grep genome size
    # while not re.search(r'^initializeGenomicElement\(g1', l):
    while not l.startswith('initializeGenomicElement(g1'):
        l = slimout.readline()
        line_number += 1
    # print(l)
    genomeSize = int(re.findall(r'[0-9]+', l)[-1])

    # same procedure for an arbitrary number of samples:
    pop_index = 0
    found_t0 = False
    while l != '':
        # go to output line "#OUT"
        # while not re.match(r'^#OUT', l) and l != '':
        while not l.startswith('#OUT') and l != '':
            l = slimout.readline()
            line_number += 1
        line = l.strip().split(' ')
        # tracking alleles:
        ## sample tracking output[annotation]:
        # #OUT:[0] 17101[1cycle] 17101[2generation] T[3tracking] p1[4subpop] 72820[5mutation id] m2[6mut type] 10000[7position] 0.02[8s] 0.5[9h] p1[10origin pop] 17099[11t0] 2[12#copy]
        if len(line) < 4:
            # print(line)
            continue
        elif line[3] == 'T' and line[4] == 'p1' and line[6] == 'm2':
            # print(line)
            assert neutral_flag == False
            # in the case of standing variation:last_id == line[5] and
            if not found_t0:
                print(" ".join(line))
                # this has to be the first tracked record, i.e. generation time is the t0 (when selection starts)
                assert int(line[1]) > int(line[11])
                # gen = 0
                # convert t0 into the absolute gen time when selection start
                t0 = int(line[1])
                gen.append(t0)
                # counts[gen] = int(line[12])
                counts.append(int(line[12]))
                found_t0 = True
                target_pos = int(line[7])
                print(f'First record of target_pos: {target_pos}. found_t0={found_t0}')
            elif last_id == line[5] and found_t0:
                try:
                    assert (target_pos == int(line[7])) or (target_pos is None), \
                        f'target_pos = {target_pos}; int(line[7]) = {int(line[7])}\n{line}'
                except Exception as e:
                    # it could be the start of another iteration
                    print(e)
                    # restart everything
                    t0 = int(line[1])
                    gen = [t0]
                    counts = [int(line[12])]
                    found_t0 = True
                    target_pos = int(line[7])
                    last_id = line[5]
                    continue
                t = int(line[1])
                gen.append(t)
                counts.append(int(line[12]))
            else:
                # then last_id != line[5]
                # restart everything
                t0 = int(line[1])
                gen = [t0]
                counts = [int(line[12])]
                found_t0 = True
                target_pos = int(line[7])
            last_id = line[5]

        # line: #OUT: 200000[1cycle] 200000[2generation] SM p1 40
        elif line[3] == 'SM' and line[4] == 'p1':
            # print(pop_index, l)
            Pop.append({})  # dict of site:count
            # record sampling time & sample size
            times.append(int(line[1]))
            sampleSize = int(line[-1])
            Sizes.append(sampleSize)

            # next next line
            l = slimout.readline()  # //
            line_number += 1
            l = slimout.readline()
            line_number += 1
            numSites.append(int(l.strip().split(' ')[-1]))
            # record positions
            l = slimout.readline().strip().split(' ')[1:]
            line_number += 1
            positions = [genomeSize * float(i) for i in l]
            # record dup:
            for pos in positions:
                if positions.count(pos) > 1:
                    dup.append(pos)

            all_positions = all_positions + positions

            Pop[pop_index] = dict.fromkeys(positions, 0)
            # start recording each sample
            for i in range(sampleSize):
                l = slimout.readline()
                line_number += 1
                ind = l.strip()
                assert len(ind) == numSites[pop_index], f'line {line_number}: unmatching length to samples' # sanity check
                for j in range(numSites[pop_index]):
                    Pop[pop_index][positions[j]] += int(ind[j])

            pop_index += 1

        elif line[3] == 'SM' and line[4] == 'p2':  # ref pop
            assert line[-1] == '1'
            # next next line
            l = slimout.readline()  # //
            line_number += 1
            l = slimout.readline()
            line_number += 1
            numSites.append(int(l.strip().split(' ')[-1]))
            # record positions
            l = slimout.readline()
            line_number += 1
            positions = l.strip().split(' ')[1:]
            positions = [genomeSize * float(i) for i in positions]

            RefPop = dict.fromkeys(positions, 1)

            # record dup:
            for pos in positions:
                if positions.count(pos) > 1:
                    dup.append(pos)

            all_positions = all_positions + positions
        # no need to read the next line bc it's all "1"s

        # read next line
        l = slimout.readline()
        line_number += 1
        if l.startswith('Mutation lost'):
            print(l)
            # just in case
            if re.search(r'[R|r]estart(ing)?', l):
                # clean slate
                restart_counter += 1
                pop_index = 0
                Pop = []  # list of dictionaries
                Sizes = []  # sample sizes
                times = []
                numSites = []  # for sanity check
                all_positions = []
                dup = []

                # for tracking allele freq
                found_t0 = False
                target_pos = None
                t0 = 0
                gen = []
                counts = []
    # eof. close file
    slimout.close()
    print('Simulation restarted', restart_counter, 'times.')

    # trim the poplist:
    all_positions = list(set(all_positions))
    all_positions.sort()  # ascending order

    # remove dups
    dup = list(set(dup))
    for dupsite in dup:
        if (dupsite == target_pos) or (dupsite == target_pos - 1):
            print(f'Target position {dupsite} is duplicated and will be removed from dataset.')
        all_positions.remove(dupsite)

    return Pop, RefPop, Sizes, times, all_positions, t0, gen, counts, target_pos


# polarize the sequence while writing output
def writeOutput(Pop, RefPop, Sizes, times, all_positions, outname, Lambda, t0, gen, counts, target_pos):
    # create file and write header:
    outfile = open(outname, 'w')
    numGen = len(counts)
    if target_pos is not None:
        outfile.write(f'##SLiM sims scaled with Lambda={Lambda}. Pop size is 10000.\n'
                      f'##Selection started at generation {t0} on position {target_pos} at'
                      f' frequency {counts[0] / 2e4} and tracked till {gen[-1]} gen in simulation.\n')  # sim.gen.scaled\tgen.ago\tcount.scaled\n
        outfile.write('##Traj.gen.ago: ' + ', '.join([str((gen[-1] - gen_i) * Lambda) for gen_i in gen]) + '\n')
        outfile.write('##Traj.count: ' + ', '.join([str(counts[i]) for i in range(numGen)]) + '\n')
    else:
        outfile.write(f'##SLiM sims scaled with Lambda={Lambda}. Pop size is 10000.\n##No Selection\n')

    # convert times to generations ago:
    times = sorted(list(set(times)))
    T = len(times)
    mostRecent = times[-1]
    times = [str((mostRecent - t) * Lambda) for t in times]
    outfile.write('##SampTimes.gen.ago: ' + ', '.join(times) + '\n')
    outfile.write('locus\t' + '\t'.join(['x%d\tn%d' % (i + 1, i + 1) for i in range(T)]) + '\n')

    # go through all the positions
    for pos in all_positions:
        if pos not in RefPop:
            # then all 1's are considered derived, print as usual
            for i in range(len(times)):
                if pos not in Pop[i]:
                    Pop[i][pos] = 0
            # make sure target_pos is written out regardless
            if pos == target_pos:
                # output line
                outfile.write('%d\t%s\n' % (pos, '\t'.join(['%d\t%d' % (Pop[j][pos], Sizes[j]) for j in range(T)])))
            elif pos == (target_pos - 1):
                # output line
                outfile.write('%d\t%s\n' % (pos, '\t'.join(['%d\t%d' % (Pop[j][pos], Sizes[j]) for j in range(T)])))
            # remove if fixed
            elif sum([Pop[j][pos] for j in range(T)]) == sum(Sizes) or sum([Pop[j][pos] for j in range(T)]) == 0:
                continue
            else:
                # output line
                outfile.write('%d\t%s\n' % (pos, '\t'.join(['%d\t%d' % (Pop[j][pos], Sizes[j]) for j in range(T)])))
        else:
            # assert Pop[ref_index][pos] == 1 #Pop no longer has ref pop
            # then all 1's are considered ancestral, flip
            for i in range(T):
                if pos not in Pop[i]:
                    Pop[i][pos] = 0
            if pos == target_pos:
                # output line
                outfile.write('%d\t%s\n' % (pos, '\t'.join(['%d\t%d' % ((Sizes[j] - Pop[j][pos]), Sizes[j]) for j in range(T)])))
            elif pos == (target_pos - 1):
                # output line
                outfile.write('%d\t%s\n' % (pos, '\t'.join(['%d\t%d' % ((Sizes[j] - Pop[j][pos]), Sizes[j]) for j in range(T)])))
            # remove if fixed
            if sum([Pop[j][pos] for j in range(T)]) == 0 or sum([Pop[j][pos] for j in range(T)]) == sum(Sizes):
                continue
            else:
                # output line
                outfile.write('%d\t%s\n' % (pos, '\t'.join(['%d\t%d' % ((Sizes[j] - Pop[j][pos]), Sizes[j]) for j in range(T)])))
    # close file
    outfile.close()


def main():
    infile = sys.argv[1]
    outfile = sys.argv[2]
    # ref_index = int(sys.argv[3])
    Lambda = int(sys.argv[3])
    # if len(sys.argv) > 4:
    #     start_freq = float(sys.argv[4])
    # else:
    #     start_freq = 20 / 2e4

    print('%s. Reading input...' % (time.ctime()))
    last = time.time()
    Pop, RefPop, Sizes, times, all_positions, t0, gen, counts, target_pos = readInput(infile)
    if (target_pos not in all_positions) and (target_pos - 1 not in all_positions):
        print(f'{target_pos} removed from all_positions.')
        # sys.exit(2)
    print('Done in %s.\nWriting output... ' % (time.time() - last))
    last = time.time()
    # writeTraj(t0, counts, trajfile, Lambda)
    writeOutput(Pop, RefPop, Sizes, times, all_positions, outfile, Lambda, t0, gen, counts, target_pos)
    print('Done in %s.' % (time.time() - last))


if __name__ == '__main__':
    main()
