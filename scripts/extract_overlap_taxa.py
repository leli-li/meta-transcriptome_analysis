#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 11 16:28:34 2025

@author: lily
"""

from pathlib import Path

# Expected structure: find_overlap/{all_find_overlap.py, kk2/, dm/, overlap/}
base_path = Path(__file__).resolve().parent
kk2_path = base_path / "kk2"
dm_path = base_path / "dm"
out_path = base_path / "overlap"
out_path.mkdir(parents=True, exist_ok=True)


def find_overlap(kk2_f, dm_f, out_f):
    kk_list = []
    dm_list = []
    overlap_list = []

    with (kk2_path / kk2_f).open("r", encoding="utf-8") as kk2:
        for line in kk2:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            taxa = line.split("\t", 1)[0]
            taxa = taxa.strip('"').replace("_", " ")
            kk_list.append(taxa)

    with (dm_path / dm_f).open("r", encoding="utf-8") as dm:
        for line in dm:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            taxa = line.split("\t", 1)[0]
            taxa = taxa.strip('"').replace("_", " ")
            dm_list.append(taxa)

    # Duplicate taxa would be written repeatedly and distort downstream analysis.
    if len(kk_list) != len(set(kk_list)):
        raise ValueError(f"Duplicate taxa found in {kk2_path / kk2_f}")
    if len(dm_list) != len(set(dm_list)):
        raise ValueError(f"Duplicate taxa found in {dm_path / dm_f}")

    kk_set = set(kk_list)
    for taxa in dm_list:
        if taxa in kk_set:
            overlap_list.append(taxa)

    with (out_path / out_f).open("w", encoding="utf-8") as out:
        for taxa in overlap_list:
            out.write(taxa + "\n")

    # 将结果写入表格
    category = Path(out_f).stem
    table_file.write(
        f"{category}\t{len(kk_list)}\t{len(dm_list)}\t{len(overlap_list)}\n"
    )

    return len(kk_list), len(dm_list), len(overlap_list)


# 创建结果表格文件并执行分析
with (out_path / "results_table.tsv").open("w", encoding="utf-8") as table_file:
    table_file.write("Category\tKK2_Count\tDM_Count\tOverlap_Count\n")

    find_overlap("bac_f.txt", "bac_f.txt", "bac_f.txt")
    find_overlap("bac_g.txt", "bac_g.txt", "bac_g.txt")
    find_overlap("bac_sp.txt", "bac_sp.txt", "bac_sp.txt")

    find_overlap("fungi_f.txt", "fungi_f.txt", "fungi_f.txt")
    find_overlap("fungi_g.txt", "fungi_g.txt", "fungi_g.txt")
    find_overlap("fungi_sp.txt", "fungi_sp.txt", "fungi_sp.txt")