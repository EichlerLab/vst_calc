import os

if not os.path.exists("./log"):
    os.makedirs("log")

localrules: all

rule all:
    input: expand("vst/{name}_{ds}_{dt}.vst.tab", name = "all_regions", ds = ["1kg", "hgdp"], dt = ["wssd", "sunk"])

rule get_vst:
    input: "data/{name}_{ds}_{dt}.genotypes.df"
    output: "vst/{name}_{ds}_{dt}.vst.tab"
    params: sge_opts = "-N get_vst -l mfree=8G", exclude_groups = "ARC D\&N CHIMPANZEE BONOBO GORILLA ORANGUTAN"
    shell:
        "python get_vst_from_table.py {input} field super_pop copy_num {output} --exclude_groups {params.exclude_groups}"
