import numpy as np
import sys
import bmws.cli




# bmws analyze results_bmws/additive_var.vcf results_bmws/additive_var.meta -d pseudohaploid -l 7 -g 1 -n 1000 >results_bmws/additive_var.results



# python bmws/bmws_change.py analyze results_bmws/additive_const.vcf results_bmws/additive_const.meta -d pseudohaploid

# this is a modified copy of the function bmws.cli.analyze_data
# the original function returns the mean selection coefficient, averaged across all generations
# the modified function averages the selection coefficient only in a specified generation interval
def analyze_data_change(args, selectionInterval):

    # original code
    meta = bmws.cli.read_meta_information(args)
    data = bmws.cli.vcf(args.vcf)
    ids = data.ids


    # modified code inserted: check that selection interval makes sense
    minSamplingGen = min(meta.values())
    maxSamplingGen = max(meta.values())
    assert (minSamplingGen <= selectionInterval[0]), (minSamplingGen, selectionInterval[0])
    assert (selectionInterval[1] <= maxSamplingGen), (selectionInterval[1], maxSamplingGen)


    # original code
    lam = 10 ** args.lam
    for snpinfo, gt in data:
        obs = bmws.cli.gt_to_obs(ids, gt, meta, args.data)
        Ne = np.full(len(obs) - 1, args.Ne)
        res, prior = bmws.cli.estimate_em(obs, Ne, lam=lam, em_iterations=args.em)


        # # original code
        # smn = np.mean(res)
        # sl1 = np.sqrt(np.mean(res * res))
        # sl2 = np.sqrt(np.mean((res - np.mean(res)) ** 2))
        # freq = np.sum(obs[:, 1]) / np.sum(obs[:, 0])

        # modified code : just take mean in interval, and also record max in and out of interval
        smn = np.mean(res[selectionInterval[0]:selectionInterval[1]])
        sl1 = -np.max(res[selectionInterval[0]:selectionInterval[1]])
        sl2 = -np.max(np.concatenate([res[:selectionInterval[0]],res[selectionInterval[1]:]]))
        freq = np.sum(obs[:, 1]) / np.sum(obs[:, 0])


        # original code
        info = [
                    str(round(freq, 3)),
                    str(round(-smn, 6)),
                    str(round(sl1, 6)),
                    str(round(sl2, 6)),
                ]

        if args.traj:
            info = info + [str(round(-f, 6)) for f in res]
            
        print(
            "\t".join(
                snpinfo
                + info
            ),
            file=args.out,
        )


def bmws_change_main (arg_list=None):

    # hack to get an additional argument into the mix for the changing seletive pressure
    selectionInterval = None
    if ('--selection_interval' in arg_list):

        # see where it is
        intIndex = arg_list.index ('--selection_interval')

        # next index should have two integers seperated by commas
        intArg = arg_list[intIndex+1]
        intList = intArg.split(',')
        assert (len(intList) == 2), len(intList)
        selectionInterval = [int(x.strip()) for x in intList]

        # nothing to see here
        del arg_list[intIndex:(intIndex+2)]
    else:
        assert (False), "Required argument missing: --selection_interval x,y"

    # parse with the defult parser
    parser = bmws.cli.get_parser()
    args = parser.parse_args(arg_list)

    if (selectionInterval is None):
        # just standard stuff
        args.func(args)
    else: 
        # analyze with the modified function
        analyze_data_change (args, selectionInterval)


def main():
    # get some args
    arg_list = sys.argv
    # remove one argument since we call it explicitly as a python script
    bmws_change_main (arg_list[1:])


if __name__ == "__main__":
    main()
