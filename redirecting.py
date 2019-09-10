from Bio import SeqIO
import os

#list_name = "RF00050.fa RF00059.fa RF00080.fa RF00162.fa RF00167.fa"
#list_name = "RF00168.fa RF00174.fa RF00234.fa RF00379.fa RF00380.fa"
#list_name = "RF00442.fa RF00504.fa RF00521.fa RF00522.fa RF00634.fa"
#list_name = "RF01054.fa RF01055.fa RF01056.fa RF01057.fa RF01482.fa"
#list_name = "RF01510.fa RF01689.fa RF01725.fa RF01727.fa RF01734.fa"
#list_name = "RF01739.fa RF01750.fa RF01767.fa RF01786.fa RF01826.fa"
#list_name = "RF01831.fa RF02680.fa RF02683.fa RF02885.fa RF02912.fa"
list_name = "RF03057.fa RF03058.fa RF03071.fa RF03072.fa"
#list_name = ['RF00379.fa','RF03072.fa']
list_name = list_name.split()

save_path1 = os.getcwd()+'/'


for fasta_name in list_name:

    save_path = save_path1 + fasta_name.split(".")[0]
    
    if not(os.path.exists(save_path)):
        os.mkdir(save_path)

    
    for seq_record in SeqIO.parse(fasta_name , "fasta"):
        part_of_name = replaceMultiple(str(seq_record.id), [".","/","-"],"_")
        file_output_name1 = fasta_name.split(".")[0] + "_" + part_of_name + '.fa'
        file_output_name = os.path.join(save_path, file_output_name1)
        file_out = open(file_output_name,'a+')
        file_out.write('> '+ str(seq_record.id) + "\n")
        file_out.write(str(seq_record.seq) + "\n")
        file_out.close()
          

