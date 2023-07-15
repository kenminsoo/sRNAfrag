from basics import *

# Description:
# Get information about annotation files.
# - Compare internal sequence with pulled sequence
# - count # of times value in column appears
# - count # of times an attribute appears

# Input gtf file, if sequence exists, this will compare the extracted one (with bedtools) to see if any errors exist in the annotation
def compare_sequence(gtf_file, output_name, sequence_feature1, sequence_feature2):
    # store data in here

    convert_df = []

    with open(gtf_file, "r") as gtf:
        for line in gtf:
            temp_line = line
            # split
            split_line = line.strip().split(sep = "\t")
            # get only sequence info
            attributes = split_line[-1].split(sep = ";")
            # strip
            attributes = list(map(str.strip, attributes))
            # combine
            attributes_split = [item.split(sep = " ") for item in attributes]
            # combine
            attributes_split = sum(attributes_split, [])
            
            # this is done because combining annotation files often results in
            # different index numbers for features

            # find index of first feature and pull sequence
            index1 = attributes_split.index(sequence_feature1)
            seq1 = attributes_split[index1 + 1].replace('"', "").upper()

            # find index of second feature and pull sequence
            index2 = attributes_split.index(sequence_feature2)
            seq2 = attributes_split[index2 + 1].replace('"', "").upper()

            # check if same length
            if len(seq1) != len(seq2):
                dist = "error"
            else:
                dist = hamming(seq1, seq2)
            temp_data = [temp_line.strip(), dist, seq1, seq2]
            convert_df.append(temp_data)
    df_for_output = pd.DataFrame(convert_df)
    df_for_output.to_csv(output_name)

# count by a certain column value
def countby_field(input_gtf, output_name, field_index):
    
    return_data = my_dictionary()
    
    with open(input_gtf, "r") as gtf, open(output_name, "w") as new:
        for line in gtf:
            split_line = line.strip().split(sep = "\t")
            value = split_line[field_index]

            if value not in return_data:
                return_data.add(value, 0)
            
            return_data[value] += 1
        
        new.write(str(return_data))

# count by an attribute
def countby_attribute(input_gtf, output_name, countby_value, skip = False, num = 1):
    
    return_data = my_dictionary()

    with open(input_gtf, "r") as gtf, open(output_name, "w") as new:
        
        iter = 1
        for line in gtf:
            if skip == True: 
                if iter <= num:
                    iter += 1
                    continue
            split_line = line.strip().split(sep = "\t")
            # get only sequence info
            attributes = split_line[-1].split(sep = ";")
            # strip
            attributes = list(map(str.strip, attributes))
            # combine
            attributes_split = [item.split(sep = " ") for item in attributes]
            # combine
            attributes_split = sum(attributes_split, [])
            # this is done because combining annotation files often results in
            # different index numbers for features
            # find index of first feature and pull sequence
            index1 = attributes_split.index(countby_value)
            type_feature = attributes_split[index1 + 1]

            if type_feature not in return_data:
                return_data.add(type_feature, 0)

            return_data[type_feature] += 1

        new.write(str(return_data))


if __name__ == '__main__':
  fire.Fire()