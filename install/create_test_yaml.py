# Will create the test yaml file

with open("install_wd.txt", "r") as wd:
    working_dir = wd.readline().strip("\n")

i = 0
print(working_dir)

# Extract the directory that sRNAfrag is installed
splitted = working_dir.split(sep = "/")

length = len(splitted)

annotation_file = "/".join(splitted[0:length - 1]) + "/scripts/H.sapiens/h.sapiens_snoRNA_final.gtf"

working = working_dir + "/test_results/working"
samples = working_dir + "/test_data"
out = working_dir + "/test_results/out"

parameters = ["      location: " + '"' + annotation_file + '"', '  working_dir: ' + '"' + working + '"', '  sample_dir: '  + '"' + samples + '"', '  out_dir: ' + '"' + out + '"']

with open("sRNA_frag_config_setup_template.yaml", "r") as config, open("sRNA_frag_config.yaml", "w") as new_config:
    for line in config:
        if '"' in line and "/" in line:
            new_config.write(parameters[i] + "\n")

            i += 1

        else:
            new_config.write(line)