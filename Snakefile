import json
import itertools

configfile: "config.yaml"
R = config["R"]

# WILDCARDS --------------------------------------------------------------------

# motifs
MOT = ["BANP", "ESR1", "GATA1","GATA2","KLF1","MYC","NR1H3", "NR1H4",
"RUNX1","RUNX2"] # "NR3C1",
# rowTypes
ROW = ["peaks"] #"motifs"
# smoother
SMT = ["none","smooth.2d_0.5", "smooth.2d_1", "smooth.2d_2"]
# differential analysis methods
DIF = glob_wildcards("code/02-dif-{x}.R").x
# data filenames
#WGT = ["weight-{}-{}-{}".format(m,r,s) for m in MOT for r in ROW for s in SMT]
#TTL = ["total-{}-{}".format(m,r) for m in MOT for r in ROW]

frg = expand("data/00-frg/frg_{mot}.rds", mot=MOT)
ttl = expand("data/01-total/total-{mot}-{row}.rds", mot=MOT, row=ROW)
wgt = expand("data/01-weight/weight-{mot}-{row}-{smt}.rds", mot=MOT, row=ROW, smt=SMT)
dif = [
    expand("outs/dif-total-{mot}-{row},{dif}.rds", mot=MOT, row=ROW, dif=DIF),
    expand("outs/dif-weight-{mot}-{row}-{smt},{dif}.rds", mot=MOT, row=ROW, smt=SMT, dif=[x for x in DIF if "chromVAR" not in x])]

# visualization
plt = []
VAL = ["dif"]
for val in VAL:
    x = glob_wildcards("code/04-plt_"+val+"-{x}.R").x
    plt += expand("plts/{val}-{plt}.pdf", val=val, plt=x)

res = {
    "frg": frg,
    "ttl": ttl,
    "wgt": wgt,
    "dif": dif
}
# SETUP ========================================================================
rule all: 
    input:
        # fragments
        #frg, 
        # counts
        ttl, wgt,
        # differential anlysis
        dif, 
        # visualization
        plt


rule get_frg:
    priority: 99
    input:  "code/00-get_frg.R",
    output: "data/00-frg/frg_{mot}.rds"
    log:    "logs/00-get_frg-{mot}.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args\
        mot={wildcards.mot} res={output}" {input[0]} {log}'''
        
rule get_ttl:
    priority: 98
    input:  "code/01-get_total.R",
            "data/00-frg/frg_{mot}.rds"
    output: "data/01-total/total-{mot}-{row}.rds"
    log:    "logs/01-total-{mot}-{row}.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args\
        mot={wildcards.mot} row={wildcards.row} res={output} frg={input[1]}" {input[0]} {log}'''

rule get_wgt:
    priority: 98
    input:  "code/01-get_weight.R",
            "data/00-frg/frg_{mot}.rds"
    output: "data/01-weight/weight-{mot}-{row}-{smt}.rds"
    log:    "logs/01-weight-{mot}-{row}-{smt}.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args\
        mot={wildcards.mot} row={wildcards.row} smt={wildcards.smt} res={output} frg={input[1]}" {input[0]} {log}'''

rule dif_ttl:
    priority: 97
    input:  "code/02-dif.R",
            "code/02-dif-{dif}.R",
            rules.get_ttl.output
    output: "outs/dif-total-{mot}-{row},{dif}.rds"
    log:    "logs/dif-total-{mot}-{row},{dif}.rds.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args\
        mot={wildcards.mot} row={wildcards.row} dif={wildcards.dif} fun={input[1]} res={output} dat={input[2]}" {input[0]} {log}'''

rule dif_wgt:
    priority: 97
    input:  "code/02-dif.R",
            "code/02-dif-{dif}.R",
            rules.get_wgt.output
    output: "outs/dif-weight-{mot}-{row}-{smt},{dif}.rds"
    log:    "logs/dif-weight-{mot}-{row}-{smt},{dif}.rds.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args\
        mot={wildcards.mot} row={wildcards.row} smt={wildcards.smt} fun={input[1]} dif={wildcards.dif} res={output} dat={input[2]}" {input[0]} {log}''' 



# VISUALIZATION ========================================================
for val in VAL:
    rule:
        priority: 49
        input:  expand("code/04-plt_{val}-{{plt}}.R", val = val), x = res[val]
        params: lambda wc, input: ";".join(input.x)
        output: expand("plts/{val}-{{plt}}.pdf", val = val)
        log:    expand("logs/plt_{val}-{{plt}}.Rout", val = val)
        shell:  '''
            {R} CMD BATCH --no-restore --no-save "--args\
            {params} {output[0]}" {input[0]} {log}'''

        
        
