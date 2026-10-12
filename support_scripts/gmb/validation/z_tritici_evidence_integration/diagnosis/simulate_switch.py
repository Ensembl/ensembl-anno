#!/usr/bin/env python3
"""Simulate replacing a backbone primary by one competitor under candidate rules (dev only)."""
import sys
import pandas as pd
df = pd.read_csv(sys.argv[1], sep="\t")
df["c_complete"] = df.c_orf.str.contains("ATG\\|STOP", regex=True)
df["p_net"] = df.p_comp - df.p_incomp; df["c_net"] = df.c_comp - df.c_incomp
rank = {"ref_cds_exact": 3, "ref_cds_ic": 2, "ref_overlap": 1}
rules = {
    "R1 fewer SR-unsupported introns": lambda d: d.c_sr_miss < d.p_sr_miss,
    "R2 R1 + complete ORF": lambda d: (d.c_sr_miss < d.p_sr_miss) & d.c_complete,
    "R3 R2 + protein net >= backbone": lambda d: (d.c_sr_miss < d.p_sr_miss) & d.c_complete & (d.c_net >= d.p_net),
    "R4 R2 + CDS >= 90% of backbone": lambda d: (d.c_sr_miss < d.p_sr_miss) & d.c_complete & (d.c_cds >= 0.9 * d.p_cds),
    "R5 R3 + CDS >= 90%": lambda d: (d.c_sr_miss < d.p_sr_miss) & d.c_complete & (d.c_net >= d.p_net) & (d.c_cds >= 0.9 * d.p_cds),
    "R6 R5 + both SR tracks agree": lambda d: (d.c_sr_miss < d.p_sr_miss) & d.c_complete & (d.c_net >= d.p_net) & (d.c_cds >= 0.9 * d.p_cds) & d.c_sources.str.contains("Scallop") & d.c_sources.str.contains("StringTie"),
}
for name, rule in rules.items():
    ok = df[rule(df)].copy()
    # one replacement per primary: fewest unsupported introns, best protein net, longest CDS
    ok = ok.sort_values(["c_sr_miss", "c_net", "c_cds"], ascending=[True, False, False]).drop_duplicates(["w", "p_rep"])
    pl = ok.p_label.map(rank).fillna(0); cl = ok.c_label.map(rank).fillna(0)
    print(f"{name:36s} switches {len(ok):4d} | exact gained {int((ok.c_exact & ~ok.p_exact).sum()):3d} "
          f"lost {int((ok.p_exact & ~ok.c_exact).sum()):3d} | label up {int((cl > pl).sum()):3d} "
          f"down {int((cl < pl).sum()):3d} same {int((cl == pl).sum()):3d} | competitor role "
          f"{ok.c_role.value_counts().to_dict()}")
