with open("At10.gtf", 'r') as file_in:
    with open("At10_new.gtf", 'w') as file_out:
        for line in file_in:
            line = line.strip('\n')
            if len(line.split('\t')) == 9:  # make sure it's data line and not header
                scaffold, source, type, start, end, score, strand, phase, attributes = line.split('\t')
                if int(start) <= int(end):
                    file_out.writelines(line + "\n")
