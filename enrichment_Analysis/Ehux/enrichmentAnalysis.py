import os, fnmatch, subprocess

modules = []
listOfFiles = os.listdir('.')
pattern = "*.ids"
for entry in listOfFiles:
        if fnmatch.fnmatch(entry, pattern):
                modules.append(entry)
modules.remove("all.ids")
print(len(modules), "modules")
print(modules)

for module in modules:
        command = "/w2/xyzhang/src/trinity/Analysis/DifferentialExpression/run_GOseq.pl --genes_single_factor " + module + " --GO_assignments ehux_go.txt --lengths ehux.seq_lens --background all.ids"
        print command
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
        process.wait()
        print "Return Code for", module, "module is", process.returncode, "\n"
