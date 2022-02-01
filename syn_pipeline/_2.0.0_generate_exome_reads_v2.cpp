// ./_2.0.0_generate_exome_reads_v2.cpp <donor_name> <synthetic_exome>.fasta <outdir>
// generates R1 and R1 fastq files that mimic Illumina reads for a given 
// input exome.fasta file. randomly samples these reads across the entire 
// exome using uniform distribution

#include <iostream>
#include <string>
#include <fstream>
#include <random>
#include <map>
#include <typeinfo>

//int read_length = 101;
std::string::size_type read_length = 101;
std::string adapter_seq = "AGATCGGAAGAGC";
int num_reads = 75900;
std::map<std::string, std::string> complement_dna { 
    {"A", "T"}, 
    {"T", "A"}, 
    {"C", "G"}, 
    {"G", "C"},
    {"N", "N"},
};

//returns reverse complement of input dna sequence
std::string reverse_complement(std::string dna_in){
    std::string dna_out = "";
    for (int i = dna_in.length()-1; i >= 0; i--){
        dna_out.append(complement_dna[dna_in.substr(i,1)]);
    }
    return dna_out;
}

int main(int argc, char** argv){
    if (argc!=5){
        std::cout << "ERROR: argc should be 5.\n";
        return 1;
    }
    std::string donor_name=argv[1];
    std::string exome_path_allele_0=argv[2];
    std::string exome_path_allele_1=argv[3];
    std::vector<std::string> exome_paths = {exome_path_allele_0, 
        exome_path_allele_1};

    std::string outdir = argv[4];
    std::cout << "outdir:" << outdir << "\n";
    std::string r1_out = outdir + "/" + donor_name +"_R1.fastq";
    std::string r2_out = outdir + "/" + donor_name +"_R2.fastq";
    std::cout << "outdir: " << outdir << "\n";
    std::cout << "r1_out: " << r1_out << "\n";
    std::cout << "r2_out: " << r2_out << "\n";

    std::ofstream r1_file;
    r1_file.open (r1_out);
    std::ofstream r2_file;
    r2_file.open (r2_out);

    for (int i= 0; i < 2; i++){
        std::string exome_path = exome_paths[i];
        std::ifstream exome_file (exome_path.c_str());
        std::string line;
        std::string exome_str = "";
        if (exome_file.is_open()){
            while ( getline (exome_file, line)){
                if(line[0]!='>'){
                    exome_str.append(line);
                }
                else{
                    std::cout << line << "\n";
                }
            }
            exome_file.close();
        }
        else {
            std::cout << "couldn't read exome_file\n";
            return 1;
        }

        std::string::size_type exome_length = exome_str.length();

        std::random_device rd; // obtain a random number from hardware
        std::mt19937 gen(rd()); // seed the generator
        std::uniform_int_distribution<std::string::size_type> distrib(0, exome_length-read_length-1); // define the range

        for (int read_num=0; read_num<num_reads; read_num++){
            // std::cout << "starting read loop\n";
            std::string::size_type start_point = distrib(gen);
            //std::cout << "start point it =  " << start_point << "\n";
            std::string read_seq = exome_str.substr((std::string::size_type)start_point, read_length);
            //std::cout << "substr done\n";

            std::string fastq_r1_l1 = "@read_";
            fastq_r1_l1.append(std::to_string(read_num));
            std::string fastq_r1_l2 = read_seq + adapter_seq;
            std::string fastq_r1_l3 = "+";
            std::string fastq_r1_l4 = std::string(fastq_r1_l2.length(), 'I');

            //std::cout << "about to write to r1\n";
            r1_file << fastq_r1_l1 << "\n" << fastq_r1_l2 << "\n" ;
            r1_file << fastq_r1_l3 << "\n" << fastq_r1_l4 << "\n" ;
            
            std::string fastq_r2_l1 = "@read_";
            fastq_r2_l1.append(std::to_string(read_num));
            std::string fastq_r2_l2 = reverse_complement(fastq_r1_l2);
            std::string fastq_r2_l3 = "+";
            std::string fastq_r2_l4 = std::string(fastq_r2_l2.length(), 'I');
            r2_file << fastq_r2_l1 << "\n" << fastq_r2_l2 << "\n" ;
            r2_file << fastq_r2_l3 << "\n" << fastq_r2_l4 << "\n" ;
        }
    }
    std::cout << "here\n";
    r1_file.close();
    r2_file.close();

    return 0;
}

