import argparse
import pandas as pd
import sys
import numpy as np

def calc_vst(group1, group2, size_thresh=30, na_value=float('nan')):
    #Calculate vst between two populations
    group1 = group1[group1 != -1]
    group2 = group2[group2 != -1]
    n1 = len(group1)
    n2 = len(group2)

    if n1 == 0 or n2 == 0 or (n1 + n2 < size_thresh):
         return na_value

    n_total = n1 + n2
    v_total = np.var(np.r_[group1, group2])
    if v_total == 0: return 0.

    v1 = np.var(group1)
    v2 = np.var(group2)
    vst = (v_total - (((v1*n1)+(v2*n2))/n_total))/v_total
    return max(vst, 0.)

def calc_region_vsts(region_df, group_col, group_combos, cp_col, size_thresh=30, na_value=float("nan")):
    region_vsts = []
    for (a, b) in group_combos:
                groupa = region_df.ix[region_df[group_col] == a, cp_col]
                groupb = region_df.ix[region_df[group_col] == b, cp_col]
                region_vsts.append(calc_vst(groupa, groupb, size_thresh=size_thresh, na_value=na_value))
    return region_vsts

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description= "Calculate Vst statistics for specified groups in a dataframe")
    parser.add_argument("long_table", help="Input dataframe in 'melted' format")
    parser.add_argument("region_column", help="Column that uniquely identifies regions")
    parser.add_argument("group_column", help="Column to use for grouping copy numbers")
    parser.add_argument("cp_column", help="Column with copy number values")
    parser.add_argument("outfile", help="Path to output dataframe")
    parser.add_argument("--size_thresh", default=30, type=int, help="Minimum total group size for vst calculation (Default: %(default)s)")
    parser.add_argument("--input_na_value", default='.', help="NA value for input data (Default: %(default)s)")
    parser.add_argument("--output_na_value", default=float('nan'), help="Default value for vst calculations that don't meet pop size threshold (Default: %(default)s)")
    parser.add_argument("--contig_col", default="chr", help="Column with chromosome info (Default: %(default)s)")
    parser.add_argument("--sex_col", default="sex", help="Column with sex info (required if 'chrX' specified in contigs, default: %(default)s)")
    parser.add_argument("--contig", help="Chromosome to get vst for (Default: all)")
    parser.add_argument("--exclude_groups", nargs="+", help="Names of groups to exclude")

    args = parser.parse_args()

    if args.input_na_value == "-1":
        input_na_value = int(args.input_na_value)
    else:
        input_na_value = args.input_na_value

    df = pd.read_table(args.long_table)

    if args.contig in ["X", "chrX"]:
        df = df.loc[df[args.sex_col].isin(["F", "f", "female"])]

    groups = df[args.group_column].unique().tolist()

    if args.exclude_groups is not None:
        groups = [group for group in groups if group not in args.exclude_groups]

    group_combos = [(a, b) for i, a in enumerate(groups) for j, b in enumerate(groups) if i < j]

    group_combo_list = ["%s_%s" % (a, b) for (a, b) in group_combos]

    regions = df[args.region_column].unique().tolist()
    df_out = pd.DataFrame(index = range(len(regions)), columns = ["region", "max_vst", "mean_vst"] + group_combo_list)

    grouped = df.groupby(args.region_column)

    for i, (region, region_df) in enumerate(grouped):
        region_vsts = calc_region_vsts(region_df, args.group_column, group_combos, args.cp_column, size_thresh=30, na_value=float("nan"))
        df_out.loc[i] = [region, max(region_vsts), np.mean(region_vsts)] + region_vsts

    df_out["chr"] = df_out["region"].map(lambda x: x.split("_")[0])
    df_out["start"] = df_out["region"].map(lambda x: x.split("_")[1])
    df_out["end"] = df_out["region"].map(lambda x: x.split("_")[2])
    df_out["name"] = df_out["region"].map(lambda x: "_".join(x.split("_")[3:]))

    df_out = df_out[["chr", "start", "end", "name", "max_vst", "mean_vst"] + group_combo_list]

    df_out = df_out.convert_objects(convert_numeric = True)

    df_out.to_csv(args.outfile, sep="\t", index=False, header=True, float_format='%.6f')
