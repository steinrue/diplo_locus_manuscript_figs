"""Pipeline to "burn-in" with msprime, pass on to SLiM, run SLiM, and parse output
Usage:
python %prog <slim_input> <condition> <from> <to> <slim_args>
"""
import sys, os, subprocess, re, time
import numpy as np
import pandas as pd
from datetime import datetime, timedelta
import msprime, pyslim, tskit

wdir = '/gpfs/data/steinruecken-lab/XiaohengStuff/Diffusion_spect/Simulations/singleLocus/'
# wdir = os.getcwd()
slimpath = '/gpfs/data/steinruecken-lab/XiaohengStuff/build/slim'


def get_seeds(seedfile, num_reps):
    Seeds = []
    if not os.path.exists(seedfile):
        seed_pool = np.empty((num_reps, 2))
        seed_pool[:,0] = np.arange(num_reps)
        seed_pool[:,1] = np.random.randint(0, 2**32 - 1, num_reps)
        np.savetxt(seedfile, seed_pool, fmt="%d", delimiter="\t", header="rep\tseed", comments="")
        enough = True
    else:
        enough = False

    with open(seedfile, 'r') as seeds:
        for l in seeds:
            if l.startswith('rep'):
                continue
            rep, seed = map(int, l.strip().split("\t"))
            assert rep == len(Seeds), f'rep=={rep}, len(Seeds)=={len(Seeds)}, l: {l}\n{Seeds}'
            Seeds.append(seed)
    seeds.close()
    enough |= (len(Seeds) >= num_reps)

    if not enough:
        extras = np.random.randint(0, 2**32 - 1, num_reps - len(Seeds))
        OG_len = len(Seeds)
        with open(seedfile, 'a') as seeds:
            for i, seed in enumerate(extras):
                seeds.write(f'{OG_len + i}\t{seed}\n')
                Seeds.append(seed)
        seeds.close()
    return Seeds


def check_slim_dirs(condition):
    # root_folder = f'{wdir}500kb_{condition}_ConstN/'
    root_folder = os.path.join(wdir, f'200kb_HC_{condition}_ConstN_t9x500gen')
    if not (os.path.exists(root_folder) and os.path.isdir(root_folder)):
        os.makedirs(root_folder)
    else:
        assert os.path.isdir(root_folder)

    for subdir in ('msout', 'trees', 'count', 'likelihood'):
        subdir_path = os.path.join(root_folder, subdir)
        if not (os.path.exists(subdir_path) and os.path.isdir(subdir_path)):
            os.makedirs(subdir_path)
        else:
            assert os.path.isdir(subdir_path)


def extract_param_from_new_eidos(eidosFile):
    """Read the eidos file to extract chr length, pop size, and mut/rec rates.
    The regex patterns only works for the new eidos scripts where constants are predefined.
    The sequence must be of type "g1", and the population (for which size will be extracted) must be named "p1"
    """
    # initialize regex
    lambda_regex = re.compile(r'defineConstant\s*\(\"LMD\"\s*,\s*(\d+)\s*\)')
    mut_regex = re.compile(r'defineConstant\s*\(\"MU\"\s*, \s*(\d+\.*\d*e*-*\d+)\s*\*\s*LMD\s*\)')
    rec_regex = re.compile(r'defineConstant\s*\(\"RHO\"\s*, \s*(\d+\.*\d*e*-*\d+)\s*\*\s*LMD\s*\)')
    pop_regex = re.compile(r'defineConstant\s*\(\"NE\"\s*, .*(\d+\.*\d*e*-*\d+)\s*/\s*LMD\s*\)')
    len_regex = re.compile(r'defineConstant\s*\(\"SEQ_LEN\"\s*, \s*(\d+\.*\d*e*-*\d+)\s*\)')
    # initialize variables
    lmd = None
    mut, rho, seq_len, N = 0, 0, 0, 0
    with open(eidosFile, 'r') as eidos:
        for l in eidos:
            # print(mut_regex.findall(l), rec_regex.findall(l), pop_regex.findall(l), len_regex.findall(l))
            if lambda_regex.search(l):
                lmd = lambda_regex.findall(l)
                assert len(lmd) == 1, f'Error in extracting scaling factor lambda; LMD={lmd};\nOG line: \"{l}\".'
                lmd = int(lmd[0])
            if mut_regex.search(l):
                # print(l)
                assert lmd is not None, f'Error: LMD should\'ve been defined before MU.'
                mut = mut_regex.findall(l)
                assert len(mut) == 1, f'Error in extracting mutation rate; mu={mut};\nOG line: \"{l}\".'
                mut = float(mut[0])
            elif rec_regex.search(l):
                # print(l)
                rho = rec_regex.findall(l)
                assert len(rho) == 1, f'Error in parsing recombination rate; rec={rho}\nOG line: {l}.'
                rho = float(rho[0])
            elif seq_len == 0 and len_regex.search(l):
                # print(l)
                seq = len_regex.findall(l)
                assert len(seq) == 1, f'Error in parsing sequence length; SEQ_LEN={seq}\nOG line: {l}.'
                seq_len = int(float(seq[0]))
            elif pop_regex.search(l):
                # print(l)
                N = pop_regex.findall(l)
                assert len(N) == 1, f'Error in parsing pop size; N={N}\nOG line: {l}.'
                N = int(float(N[0]))
            # if all vars are extracted, quit reading
            elif mut != 0 and rho != 0 and seq_len != 0 and N != 0 and lmd is not None:
                # print(f'Quit reading at {l}')
                break
    # close file in case
    eidos.close()
    # print(mut, rho, seq_len, N)
    return lmd, mut, rho, seq_len, N


# call msprime to generate a pool of 2Ne lineages of sequences of length L & mutation rate mu
# use the human-chimp demography: assume 25yr/gen; human (p2) split with chimp (p1) 5MYA --> 2e5 generations ago
def initializeHC_simple_msprime(Ne: int or float, L: int or float, mu: float, rho: float, master_seed: int):
    # initialize rates
    mut_map = msprime.RateMap(position=[0, int(L)], rate=[mu])
    rec_map = msprime.RateMap(position=[0, int(L)], rate=[rho])
    # initialize population size (Ne diploids, two pops split 2e5 generations ago)
    demog = msprime.Demography()
    demog.add_population(initial_size=int(Ne))  # name="0",
    demog.add_population(initial_size=int(Ne))  # name="1",
    demog.add_population(initial_size=int(Ne))  # name="2",
    demog.add_population_split(time=int(2e5), derived=[1, 2], ancestral=0)
    # spawn 2 seeds from the given seed
    rng = np.random.RandomState(master_seed)
    anc_seed, mut_seed = rng.randint(0, 2**32 - 1, 2)
    # simulate trees
    mstrees = msprime.sim_ancestry(samples={1: int(Ne), 2: int(Ne)},
                                   demography=demog, random_seed=anc_seed, recombination_rate=rec_map)
    # use pyslim to annotate lineages
    mstrees = pyslim.annotate(mstrees, model_type="WF", tick=int(2e5-1), stage="late")
    # add neutral mutations
    m1_model = msprime.SLiMMutationModel(type=1)
    mstrees = msprime.sim_mutations(mstrees, rate=mut_map, model=m1_model,
                                    keep=True, random_seed=mut_seed)

    # add metadata
    temp_metadata = mstrees.tables.metadata
    temp_metadata["SLiM"]["model_type"] = "WF"
    mstrees.tables.metadata = temp_metadata
    # done
    return mstrees


# frac_regex = re.compile(r'(\.\d+)')
s_regex = re.compile(r's(\.\d+)')
h_regex = re.compile(r'h(\.{0,1}\d+)')
f0_regex = re.compile(r'stdVar(\.\d+)')


def get_slim_command(slim_input, condition, rep, seed, tree_file, sample_size): #, s:float=0., h: float=0.5 <-- can be extracted from cond
    global wdir, slimpath
    # tree_file = f'200kb_HC_{COND}_ConstN_t9x500gen/trees/rep{rep}_seed{rep_seed}.trees'
    slim_output = f'{wdir}200kb_HC_{condition}_ConstN_t9x500gen/msout/{condition}_ConstN_t9x500gen_n40_rep{rep}.msout'

    if 'neut' in condition.lower():
        slim_command = f'{slimpath} -t -s {seed} -d input_mstrees=\\\"{tree_file}\\\" -d sample_size={sample_size} {slim_vars} {slim_input} > {slim_output}'
    # selection
    else:
        # extract parameters
        S = s_regex.findall(condition)
        H = h_regex.findall(condition)
        # sanity check
        assert len(S) > 0 and len(H) > 0, f'Cannot parse coefficients from condition \"{condition}\".'
        S = float(S[0])
        H = float(H[0])
        if 'sdtVar-fromSFS' in condition:
            print(condition, S, H)
            assert 'stdVarSFS' in slim_input, f'Not input for standing variation: {slim_input}.'
            slim_command = f'{slimpath} -t -s {seed} -d input_mstrees=\"\'{tree_file}\'\" -d sample_size={sample_size} -d var_s={S} -d var_h={H} {slim_input} > {slim_output}'
        elif 'stdVar' in condition and 'SFS' not in condition:
            assert 'stdVar' in slim_input, f'Not input for standing variation: {slim_input}.'
            f0 = f0_regex.findall(condition)
            assert len(f0) > 0, f'Cannot extract f_0 MAF threshold from condition \"{condition}\".'
            f0 = float(f0[0])
            assert 0 < f0 < 1, f'Initial frequency threshold must be between 0 and 1: f0={f0}.'
            slim_command = f'{slimpath} -t -s {seed} -d input_mstrees=\"\'{tree_file}\'\" -d sample_size={sample_size} -d var_s={S} -d var_h={H} -d thres_freq={f0} {slim_input} > {slim_output}'
        else:
            slim_command = f'{slimpath} -t -s {seed} -d input_mstrees=\"\'{tree_file}\'\" -d sample_size={sample_size} -d var_s={S} -d var_h={H} {slim_input} > {slim_output}'
    return slim_command, slim_output


def print_time(duration):
    return time.strftime("%H hr %M min %S sec", time.gmtime(duration))


def print_timedelta(duration: timedelta):
    output_str = f'{duration.days} days {print_time(duration.total_seconds())}'
    return output_str


def track_selPos_in_msout(filename):
    # initialize
    Samples = []  # list of dictionaries
    Sizes = []  # sample sizes
    times = []  # sample times
    numSites = []  # for sanity check
    all_positions = []
    dup = []
    # for tracking allele freq
    last_id = '0'
    t0 = 0
    gen = []
    counts = []

    restart_counter = 0
    target_pos = None
    line_number = 0

    # open file
    slimout = open(filename, 'r')
    l = slimout.readline(); line_number += 1

    # grep genome size
    # while not re.search(r'^initializeGenomicElement\(g1', l):
    while not l.startswith('initializeGenomicElement(g1'):
        l = slimout.readline(); line_number += 1
    # print(l)
    genomeSize = int(re.findall(r'[0-9]+', l)[-1])

    # same procedure for an arbitrary number of samples:
    pop_index = 0
    found_t0 = False
    T_line = ''
    while l != '':
        # go to output line "#OUT"
        # while not re.match(r'^#OUT', l) and l != '':
        while not l.startswith('#OUT') and (l != ''):
            l = slimout.readline(); line_number += 1
        if l == '':
            break
        # now l == "#OUT: 199999 199999 T p1 740 m1 50045 0 0.5 p-1 -2 6"
        # or l == "#OUT: [gen] T"
        line = l.strip().split(' ')
        # tracking alleles:
        ## sample tracking output[annotation]:
        # #OUT:[0] 17101[1cycle] 17101[2generation] T[3tracking] p1[4subpop] 72820[5mutation id] m2[6mut type] 10000[7position] 0.02[8s] 0.5[9h] p1[10origin pop] 17099[11t0] 2[12#copy]
        if len(line) < 4: # <-- most likely eof
            print(line)
            continue
        elif line[3] == 'T' and line[4] == 'p1' and line[6] == 'm2':
            T_line = [str(line_number)+":"] + line
            # in the case of standing variation:
            if not found_t0:  # last_id == line[5] and
                # print(f'found_t0 == {found_t0}, target_pos == {target_pos}')
                # this has to be the first tracked record, i.e. generation time is the t0 (when selection starts)
                assert int(line[1]) > int(line[11])
                assert len(gen) == 0 and len(counts) == 0, f'len(gen)={len(gen)}, len(counts)={len(counts)}'
                # convert t0 into the absolute gen time when selection start
                target_pos = int(line[7])
                t0 = int(line[1])
                gen.append(t0)
                counts.append(int(line[12]))
                found_t0 = True
                print(f'First record of target_pos: {target_pos}. found_t0={found_t0}')
            elif last_id == line[5] and found_t0:
                # print(f'last_id == line[5] == {last_id} and found_t0 == {found_t0}')
                try:
                    assert (target_pos == int(line[7])) or (target_pos is None), \
                        f'target_pos = {target_pos}; int(line[7]) = {int(line[7])}\nline {line_number}: {" ".join(line)}'
                except Exception as e:
                    # it could be the start of another iteration
                    print(line_number, e)
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
                target_pos = int(line[7])
                # print(f'last_id == line[5] == {last_id} and found_t0 == {found_t0}')
            elif found_t0:
                # then last_id != line[5]
                assert last_id != line[5], f'last_id={last_id} line[5]={line[5]}'
                # restart everything
                t0 = int(line[1])
                gen = [t0]
                counts = [int(line[12])]
                found_t0 = True
                target_pos = int(line[7])
            else:
                print("fell through...")
                print(T_line)
                print(target_pos, found_t0, last_id, t0, gen, counts)
            last_id = line[5]
            print(f'target_pos: {target_pos}. found_t0={found_t0}')
        # line: #OUT: 200000[1cycle] 200000[2generation] SM p1 40
        elif line[3] == 'SM' and line[4] == 'p1':
            # print(pop_index, l)
            Samples.append(0)  # list of count
            # sanity check:
            assert len(Samples) == (pop_index + 1), f'len(Samples)={len(Samples)}; pop_index={pop_index};\n\"{line}\"'
            # record sampling time & sample size
            times.append(int(line[1]))
            sampleSize = int(line[-1])
            Sizes.append(sampleSize)
            # next next line
            l = slimout.readline(); line_number += 1  # //
            l = slimout.readline(); line_number += 1
            numSites.append(int(l.strip().split(' ')[-1]))
            # record positions
            l_pos = slimout.readline().strip().split(' ')[1:]
            line_number += 1
            positions = [int(genomeSize * float(i)) for i in l_pos]
            # get idx for target pos
            print(target_pos, found_t0)
            if target_pos is not None:
                if (target_pos -1) in positions:
                    target_idx = positions.index(target_pos -1)
                elif target_pos in positions:
                    target_idx = positions.index(target_pos)
                else:
                    # not observed, skip
                    pop_index += 1
                    continue
            else:
                raise ValueError(f'target_pos should exist/be recorded before reading line {line_number}:'
                                 f'{" ".join(line)}\n{" ".join(T_line)}\n'
                                 f'{l}\n{filename}')
            # record dup:
            for pos in positions:
                if positions.count(pos) > 1:
                    dup.append(pos)

            all_positions = all_positions + positions

            # start recording *the* selected site
            for i in range(sampleSize):
                l = slimout.readline(); line_number += 1
                ind = l.strip() # single hap
                try:
                    assert len(ind) == numSites[pop_index], f'len(ind) = {len(ind)}, ' \
                                                        f'pop_index = {pop_index}'
                except Exception as e:
                    print(e)
                    print(f'len(ind) = {len(ind)}, pop_index = {pop_index}, len(numSites) = {len(numSites)}')
                    sys.exit(2)
                Samples[pop_index] += int(ind[target_idx])
            pop_index += 1

        elif line[3] == 'SM' and line[4] == 'p2':  # ref pop
            assert line[-1] == '1'
            # next next line
            l = slimout.readline(); line_number += 1  # //
            l = slimout.readline(); line_number += 1
            # numSites.append(int(l.strip().split(' ')[-1]))
            # record positions
            l = slimout.readline(); line_number += 1
            positions = l.strip().split(' ')[1:]
            positions = [int(genomeSize * float(i)) for i in positions]

            # if target site is 1 here, flip samples
            if (target_pos -1) in positions:
                Samples = [Sizes[i] - count for i, count in enumerate(Samples)]

            # record dup:
            for pos in positions:
                if positions.count(pos) > 1:
                    dup.append(pos)

            all_positions = all_positions + positions
            # no need to read the next line bc it's all "1"s

        # read next line
        l = slimout.readline(); line_number += 1

        if l.startswith('Mutation lost'):
            print(l)   # <-- sometimes it can be "Mutation lost. Terminating simulation."
            if 'estart' in l:
                # clean slate
                restart_counter += 1
                pop_index = 0
                Samples = []  # list of dictionaries
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

    dup = list(set(dup))
    # check if target_pos is dup
    if target_pos is not None:
        if target_pos in dup:
            print(f'Target position {target_pos} is duplicated and will be removed from dataset.')
            target_pos = None
        elif target_pos - 1 in dup:
            print(f'Target position {target_pos - 1} is duplicated and will be removed from dataset.')
            target_pos = None
    # remove dups
    for dupsite in dup:
        all_positions.remove(dupsite)

    return target_pos, gen, counts, all_positions, Samples, Sizes


def main():
    # seed_file = sys.argv[1]
    num_reps = 200
    slim_input = sys.argv[1]
    COND = sys.argv[2]
    from_rep, to_rep = int(sys.argv[3]), int(sys.argv[4])
    masks = sys.argv[5].split(',')
    # slim_vars = sys.argv[5]

    # retrieve parameters
    LMD, MU, RHO, LEN, NE = extract_param_from_new_eidos(slim_input)

    # make dirs
    check_slim_dirs(condition=COND)

    # read seeds
    seed_file = f'{wdir}slim4_input/{(LEN/1e3):g}kb_HC_{COND}_ConstN_t9x500gen_simSeeds.txt'
    Seeds = get_seeds(seed_file, num_reps)
    print(len(Seeds))
    # sanity check
    assert 0 <= from_rep <= to_rep <= num_reps, f'start={from_rep}, end={to_rep}, num_reps={num_reps}'

    # go through reps
    rerun_count = 0
    for rep in range(from_rep, to_rep + 1):
        rep_seed = Seeds[rep]
        # tree file name
        tree_file = f'{(LEN/1e3):g}kb_HC_{COND}_ConstN_t9x500gen/trees/rep{rep}_seed{rep_seed}.trees'
        # slim stuff
        slim_command, slim_msout = get_slim_command(slim_input, COND, rep, rep_seed, tree_file, sample_size=40)
        last = datetime.now()
        run_rep = True
        new_run = False  # <-- if True, then remove old files (bc it'd be a new sim)
        while run_rep:
            print(f'rep {rep} uses master seed {rep_seed} for both ms and SLiM.')
            last = datetime.now()
            # print(tree_file)
            ## decide whether to skip:
            if not os.path.exists(tree_file):
                # run msprime
                init_trees = initializeHC_simple_msprime(NE, LEN, MU, RHO, master_seed=rep_seed)
                # output trees
                init_trees.dump(tree_file)
                new_run |= True
            elif os.path.getsize(tree_file) == 0:
                # run msprime
                init_trees = initializeHC_simple_msprime(NE, LEN, MU, RHO, master_seed=rep_seed)
                # output trees
                init_trees.dump(tree_file)
                new_run |= True
            else:
                print(f'{tree_file} exists.')
            # run slim
            slim_command, slim_msout = get_slim_command(slim_input, COND, rep, rep_seed, tree_file, sample_size=40)
            if not os.path.exists(slim_msout) or not os.path.exists(tree_file):
                slim_stdout = subprocess.check_output(slim_command, shell=True)
                print(slim_stdout.decode())
                new_run |= True
            ## check if re-run needed:
            target_pos, gen, counts, all_positions, Samples, Sizes = track_selPos_in_msout(slim_msout)
            if (target_pos is None) or (len(Samples) < 9):
                # this could happen when sim ends prematurely (mutation lost)
                # rerun
                run_rep = True
                new_run |= True
            elif (target_pos not in all_positions) and (target_pos - 1 not in all_positions):
                # rerun
                print(f'Neither {target_pos} nor {target_pos - 1} is in all_positions. Rerun.')
                run_rep = True
                new_run |= True
            else:  # position exist and mutation not lost
                run_rep = False
                for mask in masks:
                    if mask == '2morePoly':
                        # get freq
                        freqs = np.array(Samples) / np.array(Sizes)
                        is_poly = (freqs > 0) & (freqs < 1)
                        num_poly = is_poly.sum()
                        # re-run if <=2
                        run_rep |= (num_poly <= 2)
                        new_run |= (num_poly <= 2)
                    elif "MAF" in mask:
                        minMAF = float(mask[3:])
                        poolfreq = sum(Samples) / sum(Sizes)
                        # fold
                        if poolfreq > 0.5:
                            poolfreq = 1 - poolfreq
                        run_rep |= (poolfreq <= minMAF)
                        new_run |= (poolfreq <= minMAF)
                    else:
                        raise ValueError(f'cannot recognize mask {mask}')
            print(run_rep)
            if run_rep:
                # check out this run?
                print(f'rep {rep}\t{rep_seed}\tpos {target_pos}:\t{Samples}')
                # remove old files (?)
                if os.path.exists(tree_file):
                    os.remove(tree_file)
                if os.path.exists(slim_msout):
                    os.remove(slim_msout)
                # increment
                rerun_count += 1
                rep_seed += 1
                ## tree file name
                tree_file = f'{(LEN/1e3):g}kb_HC_{COND}_ConstN_t9x500gen/trees/rep{rep}_seed{rep_seed}.trees'
                ## slim stuff
                slim_command, slim_msout = get_slim_command(slim_input, COND, rep, rep_seed, tree_file, sample_size=40)
                new_run |= True
                print(f'Re-run SLiM for rep {rep} with seed {rep_seed}')
            else:
                print(f'rep {rep}\t{rep_seed}\tpos {target_pos}:\t{Samples}')

        print(f'rep{rep} rerun pipeline {rerun_count} times\t{time.ctime()}\t'
              f'{print_timedelta(datetime.now() - last)}\n\n')

        # remove files associated with the old rep
        parsed_count = f'{(LEN/1e3):g}kb_HC_{COND}_ConstN_t9x500gen/count/HC_{COND}_ConstN_t9x500gen_n40_rep{rep}.count'

        if new_run:
            print(f'Overwrite old simulation for rep {rep}')
            if os.path.exists(parsed_count):
                os.remove(parsed_count)
            # just outsource then
            ls_command = f'ls {wdir}{(LEN/1e3):g}kb_HC_{COND}_ConstN_t9x500gen/likelihood/*_rep{rep}_* | wc -l '
            file_count = subprocess.check_output(ls_command, shell=True)
            if int(file_count) > 0:
                rm_command = f'rm {wdir}{(LEN/1e3):g}kb_HC_{COND}_ConstN_t9x500gen/likelihood/*_rep{rep}_*'
                rm_stdout = subprocess.check_output(rm_command, shell=True)
                print(rm_stdout.decode())

        # now let's parse
        if not os.path.exists(parsed_count):
            # call parsing script
            parsing_script = f'{wdir}parseSlim4_TrackSel_var-t.py'
            parsing_command = f'python {parsing_script} {slim_msout} {parsed_count} {LMD}'
            parse_stdout = subprocess.check_output(parsing_command, shell=True)
            print(parse_stdout.decode())
        elif (os.path.getsize(parsed_count) == 0) or ('reparse' in sys.argv):
            # call parsing script
            parsing_script = f'{wdir}parseSlim4_TrackSel_var-t.py'
            parsing_command = f'python {parsing_script} {slim_msout} {parsed_count} {LMD}'
            parse_stdout = subprocess.check_output(parsing_command, shell=True)
            print(parse_stdout.decode())
        else:
            print(f'{parsed_count} exists.')


if __name__ == "__main__":
    main()