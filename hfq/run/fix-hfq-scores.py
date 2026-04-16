#!/usr/bin/env python
import os
import pandas as pd

def main():
    columns = ["logit","leader zscore"]
    for txt in os.listdir("hfq-scores"):
        table = pd.read_csv(os.path.join("hfq-scores",txt),sep="\t",index_col=0)
        table.columns = columns
        table.to_csv("hfq-scores.fixed/" + txt, sep="\t")

if __name__ == "__main__":
    main()
