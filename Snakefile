import json
import itertools

configfile: "config.yaml"
R = config["R"]

############### real datasets ################
# WILDCARDS --------------------------------------------------------------------
# motifs
MOT = ["BANP", "ESR1", "GATA1","GATA2","KLF1","NR1H3", "NR1H4",
"RUNX1"] # "NR3C1", "MYC", "RUNX2"
# rowTypes
ROW = ["peaks"] #"motifs"

# smoother
#SMT = ["none","smooth.2d_0.25", "smooth.2d_0.5,", "smooth.2d_0.75", "smooth.2d_1"]
SMT = ["none", "smooth.2d_1"]
PKW = ["none", "loess", "lm", "wlm"]
# differential analysis methods
DIF = glob_wildcards("code/02-dif-{x}.R").x
PKD = glob_wildcards("code/03-pkd-{x}.R").x

# RESULTS --------------------------------------------------------------------
dat = expand("/mnt/plger/jwang/data/dat/00-frg/{mot}.rds", mot=MOT)
dat_ttl = expand("/mnt/plger/jwang/data/dat/01-total/total-{mot}-{row}.rds", mot=MOT, row=ROW)
dat_wgt = expand("/mnt/plger/jwang/data/dat/01-weight/weight-{mot}-{row}-{smt}-{pkw}.rds", mot=MOT, row=ROW, smt=SMT, pkw=PKW)
dat_dif = [
    expand("outs/dat/dif-total-{mot}-{row},{dif}.rds", mot=MOT, row=ROW, dif=DIF),
    expand("outs/dat/dif-weight-{mot}-{row}-{smt}-{pkw},{dif}.rds", mot=MOT, row=ROW, smt=SMT, pkw=PKW, dif=[x for x in DIF if "chromVAR" not in x])]

dat_pkd = [
    expand("outs/dat/pkd-total-{mot}-{row},{pkd}.rds", mot=MOT, row=ROW, pkd=PKD),
    expand("outs/dat/pkd-weight-{mot}-{row}-{smt}-{pkw},{pkd}.rds", mot=MOT, row=ROW, smt=SMT, pkw=PKW, pkd=PKD)]


############### simulated datasets ################
# WILDCARDS --------------------------------------------------------------------
MTF = ["CEBPB", "MAZ", "ZNF143"]
EFT = ["0", "0.25", "0.5", "1", "3"]

# RESULTS --------------------------------------------------------------------
sim = expand("/mnt/plger/jwang/data/sim/00-frg/{mtf},{eft}.rds", mtf=MTF, eft=EFT)
sim_ttl = expand("/mnt/plger/jwang/data/sim/01-total/total-{mtf},{eft},{row}.rds", mtf=MTF, eft=EFT, row=ROW)
sim_wgt = expand("/mnt/plger/jwang/data/sim/01-weight/weight-{mtf},{eft},{row}-{smt}-{pkw}.rds", mtf=MTF, eft=EFT, row=ROW, smt=SMT, pkw=PKW)


dat_res = {
    "dat": dat,
    "ttl": dat_ttl,
    "wgt": dat_wgt,
    "dif": dat_dif,
    "pkd": dat_pkd
}

sim_res = {
     "sim": sim,
     "ttl": sim_ttl,
     "wgt": sim_wgt
}

# visualization
plt = []
VAL = dat_res.keys()
for val in VAL:
    x = glob_wildcards("code/04-plt_"+val+"-{x}.R").x
    plt += expand("plts/dat/{val}-{plt}.pdf", val=val, plt=x)

# SETUP ========================================================================
rule all: 
    input:
        [x for x in dat_res.values()], plt,
        [x for x in sim_res.values()]

# real datasets
rule get_dat:
    priority: 99
    input:  "code/00-get_dat.R",
    output: "/mnt/plger/jwang/data/dat/00-frg/{mot}.rds"
    log:    "logs/00-get_dat-{mot}.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args\
        mot={wildcards.mot} res={output}" {input[0]} {log}'''
        
rule dat_ttl:
    priority: 98
    input:  "code/01-dat_total.R",
            "/mnt/plger/jwang/data/dat/00-frg/{mot}.rds"
    output: "/mnt/plger/jwang/data/dat/01-total/total-{mot}-{row}.rds"
    log:    "logs/01-total-{mot}-{row}.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args\
        mot={wildcards.mot} row={wildcards.row} res={output} frg={input[1]}" {input[0]} {log}'''

rule dat_wgt:
    priority: 98
    input:  "code/01-dat_weight.R",
            "/mnt/plger/jwang/data/dat/00-frg/{mot}.rds"
    output: "/mnt/plger/jwang/data/dat/01-weight/weight-{mot}-{row}-{smt}-{pkw}.rds"
    log:    "logs/dat_weight-{mot}-{row}-{smt}-{pkw}.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args\
        mot={wildcards.mot} row={wildcards.row} smt={wildcards.smt} pkw={wildcards.pkw} res={output} frg={input[1]}" {input[0]} {log}'''

rule dif_ttl:
    priority: 97
    input:  "code/02-dif.R",
            "code/02-dif-{dif}.R",
            rules.dat_ttl.output
    output: "outs/dat/dif-total-{mot}-{row},{dif}.rds"
    log:    "logs/dat_dif-total-{mot}-{row},{dif}.rds.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args\
        mot={wildcards.mot} row={wildcards.row} dif={wildcards.dif} fun={input[1]} res={output} dat={input[2]}" {input[0]} {log}'''

rule dif_wgt:
    priority: 97
    input:  "code/02-dif.R",
            "code/02-dif-{dif}.R",
            rules.dat_wgt.output
    output: "outs/dat/dif-weight-{mot}-{row}-{smt}-{pkw},{dif}.rds"
    log:    "logs/dif-weight-{mot}-{row}-{smt}-{pkw},{dif}.rds.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args\
        mot={wildcards.mot} row={wildcards.row} smt={wildcards.smt} pkw={wildcards.pkw} fun={input[1]} dif={wildcards.dif} res={output} dat={input[2]}" {input[0]} {log}''' 

rule pkd_ttl:
    priority: 96
    input:  "code/03-pkd.R",
            "code/03-pkd-{pkd}.R",
            rules.dat_ttl.output
    output: "outs/dat/pkd-total-{mot}-{row},{pkd}.rds"
    log:    "logs/dat_pkd-total-{mot}-{row},{pkd}.rds.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args\
        mot={wildcards.mot} row={wildcards.row} pkd={wildcards.pkd} fun={input[1]} res={output} dat={input[2]}" {input[0]} {log}'''

rule pkd_wgt:
    priority: 96
    input:  "code/03-pkd.R",
            "code/03-pkd-{pkd}.R",
            rules.dat_wgt.output
    output: "outs/dat/pkd-weight-{mot}-{row}-{smt}-{pkw},{pkd}.rds"
    log:    "logs/dat_pkd-weight-{mot}-{row}-{smt}-{pkw},{pkd}.rds.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args\
        mot={wildcards.mot} row={wildcards.row} smt={wildcards.smt} pkw={wildcards.pkw} fun={input[1]} pkd={wildcards.pkd} res={output} dat={input[2]}" {input[0]} {log}''' 


# simulation
rule get_sim:
    priority: 99
    input:  "code/00-get_sim.R",
    output: "/mnt/plger/jwang/data/sim/00-frg/{mtf},{eft}.rds"
    log:    "logs/00-get_sim-{mtf},{eft}.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args\
        mtf={wildcards.mtf} eft={wildcards.eft} res={output}" {input[0]} {log}'''

rule sim_ttl:
    priority: 98
    input:  "code/01-sim_total.R",
            rules.get_sim.output
    output: "/mnt/plger/jwang/data/sim/01-total/total-{mtf},{eft},{row}.rds"
    log:    "logs/01-total-{mtf},{eft},{row}.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args\
        mtf={wildcards.mtf} eft={wildcards.eft} row={wildcards.row} res={output} frg={input[1]}" {input[0]} {log}'''

rule sim_wgt:
    priority: 98
    input:  "code/01-sim_weight.R",
            rules.get_sim.output
    output: "/mnt/plger/jwang/data/sim/01-weight/weight-{mtf},{eft},{row}-{smt}-{pkw}.rds"
    log:    "logs/dat_weight-{mtf},{eft},{row}-{smt}-{pkw}.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args\
        mtf={wildcards.mtf} eft={wildcards.eft} row={wildcards.row} smt={wildcards.smt} pkw={wildcards.pkw} res={output} frg={input[1]}" {input[0]} {log}'''


# VISUALIZATION ========================================================
for val in VAL:
    rule:
        priority: 49
        input:  expand("code/04-plt_{val}-{{plt}}.R", val = val), x = dat_res[val]
        params: lambda wc, input: ";".join(input.x)
        output: expand("plts/dat/{val}-{{plt}}.pdf", val = val)
        log:    expand("logs/plt_{val}-{{plt}}.Rout", val = val)
        shell:  '''
            {R} CMD BATCH --no-restore --no-save "--args\
            {params} {output[0]}" {input[0]} {log}'''

        
        
