                        ##########################
                        ### Version 0.1        ###
                        ### Amplicon Pipe-line ###
                        ##########################

## input a file as:
## i1_S10_L001_R1_001.fastq.gz,i1_S10_L001_R2_001.fastq.gz,i1,Forward,Reverse
## i4_S13_L001_R1_001.fastq.gz,i4_S13_L001_R2_001.fastq.gz,i4,Forward,Reverse
##
## more strict filters in the mergepairs command: maxdifs=0, fastq_pctid: 100 and fastq_minimal Q>= 35
########################################################
#!/bin/bash
#Set Script Name variable
SCRIPT=`basename ${BASH_SOURCE[0]}`
#Initialize variables to default values.
CUT_F=250
CUT_R=150
MERGEID=100
CLEAN=0
TRIMBYCUT=1
PESE=1
WINQ="4:30"
MINL=75
#Set fonts for Help.
NORM=`tput sgr0`
BOLD=`tput bold`
REV=`tput smso`
#Help function
function HELP {
  echo -e \\n"Help documentation for ${BOLD}${SCRIPT}.${NORM}"\\n
  echo -e "${REV}Basic usage:${NORM} ${BOLD}$SCRIPT [FLAGS] mapping_file${NORM}"\\n
  echo "Command line switches are optional. The following switches are recognized."
  echo "${REV}-f${NORM}  --Sets the value for option ${BOLD}trimm_forward_read (int)${NORM}. Default is ${BOLD}$CUT_F${NORM}."
  echo "${REV}-r${NORM}  --Sets the value for option ${BOLD}trimm_reverse_read (int)${NORM}. Default is ${BOLD}$CUT_R${NORM}."
  echo "${REV}-i${NORM}  --Sets the value for option ${BOLD}percentageID_merge (int)${NORM}. Default is ${BOLD}$MERGEID${NORM}."
  echo "${REV}-c${NORM}  --Sets the value for option ${BOLD}clean intermediate files (Bool)${NORM}. Default is ${BOLD}$CLEAN${NORM}."
  echo "${REV}-l${NORM}  --Sets the value for option ${BOLD}initial trimming by Lenght(Bool)[byLenght=1, Quality(trimmomatic)=0]${NORM}. Default is ${BOLD}$TRIMBYCUT${NORM}."
  echo "${REV}-p${NORM}  --Sets the value for option ${BOLD}Analyse PE or SE (Bool)[PE=1 SE=0], if a PE=0 then a reference must be given as: ref_otu.fasta${NORM}. Default is ${BOLD}$PESE${NORM}."
  echo "${REV}-q${NORM}  --Sets the value for option ${BOLD}Trimm by Lenght(-l 0): WindowSize:Quality (string int:int)${NORM}. Default is ${BOLD}$WINQ${NORM}."
  echo "${REV}-m${NORM}  --Sets the value for option ${BOLD}Trimm by Lenght(-l 0): Min Size Amplicon (int)${NORM}. Default is ${BOLD}$MINL${NORM}."
  echo -e "${REV}-h${NORM}  --Displays this helpful message. No further functions are performed."\\n
  echo -e "Example: ${BOLD}$SCRIPT -f 200 -r 150 -i 98 -c 1 -l 1 mapping_file${NORM}"\\n
  echo -e "Example: ${BOLD}$SCRIPT -f 200 -r 150 -i 98 -p 0 mapping_file${NORM} *in the folder a file named: ref_otu.fasta"\\n
  exit 1
}

#Check the number of arguments. If none are passed, print help and exit.
NUMARGS=$#
echo -e \\n"Number of arguments: $NUMARGS"
if [ $NUMARGS -eq 0 ]; then
  HELP
fi

### Start getopts code ###
while getopts :f:r:i:c:l:p:m:q:h FLAG; do
  case $FLAG in
    f)
      CUT_F=$OPTARG
      echo "-cutF used: $OPTARG"
      ;;
    r)
      CUT_R=$OPTARG
      echo "-cutR used: $OPTARG"
      ;;
    i)
      MERGEID=$OPTARG
      echo "-mergeid used: $OPTARG"
      ;;
    c)
      CLEAN=$OPTARG
      echo "-clean used: $OPTARG"
      ;;
    p)
      PESE=$OPTARG
      echo "-SEorPE used: $OPTARG"
      ;;
    l)
      TRIMBYCUT=$OPTARG
      echo "-TrimbyLenght used: $OPTARG"
      ;;
    q)
      WINQ=$OPTARG
      echo "-Quality:Window used: $OPTARG"
      ;;
    m)
      MINL=$OPTARG
      echo "-Min. Length amplicon  used: $OPTARG"
      ;;
    h)
      HELP
      ;;
    \?) #unrecognized option - show help
      echo -e \\n"Option -${BOLD}$OPTARG${NORM} not allowed."
      HELP
      ;;
  esac
done
shift $((OPTIND-1))  #moving on with the arguments

echo "cutf:$CUT_F cutR:$CUT_R Mergeid:$MERGEID clean:$CLEAN"
#################################################################################################
#=Some functions#
count_fastq(){
seq=`awk '{s++}END{print s/4}' $1`
echo $seq
}
count_fasta(){
seq=`grep ">" -c $1`
echo $seq
}
#===============#
rm Stat.txt #Deleting previous summary files #####

############################################################################################################################## SE MODE ##################################
if [ $PESE -eq 0 ]
        then
        echo running in SE!!!! is it right?
        for i in $(cat $1);
                do
                array=(${i//,/ })
                java -jar /home/dich/software/Trimmomatic-0.38/trimmomatic-0.38.jar SE ${array[0]} ${array[2]}.QfiltrimSE.fastq ILLUMINACLIP:/home/dich/software/Trimmomatic-0.38/adapters/TruSeq2-SE.fa:2:30:10 LEADING:28 TRAILING:28 SLIDINGWINDOW:$WINQ MINLEN:$MINL
                cutadapt --format fastq -g ${array[3]} -o ${array[2]}.QfiltrimSE_noAdapt.fastq ${array[2]}.QfiltrimSE.fastq
                /home/dich/software/usearch_10.0.240_lin32 -fastx_uniques ${array[2]}.QfiltrimSE_noAdapt.fastq -fastaout ${array[2]}.QfiltrimSE_noAdapt_uniq.fasta -relabel Uniq
                /home/dich/software/usearch_10.0.240_lin32 -otutab ${array[2]}.QfiltrimSE_noAdapt.fastq -otus  ${array[2]}.QfiltrimSE_noAdapt_uniq.fasta -id 1 -otutabout otutab_${array[2]}.txt -mapout map_${array[2]}.txt -notmatched unmapped_${array[2]}.fa -dbmatched otus_with_sizes_${array[2]}.fa -sizeout
                /home/dich/software/usearch_10.0.240_lin32 -otutab ${array[2]}.QfiltrimSE_noAdapt.fastq -otus ref_otu.fasta -id 1 -otutabout otutab_REFF_${array[2]}.txt -mapout map_REFF_${array[2]}.txt -notmatched unmapped_REFF_${array[2]}.fa -dbmatched otus_with_sizes_REFF_${array[2]}.fa -sizeout

        done
echo SE running mode......... Done.
        if [ $CLEAN -eq 1 ]
                then
                rm *_noAdapt_uniq.fasta *noAdapt.fastq unmapped_*.fa map_*.txt
                echo Clean..... Done.
        fi
        exit
fi

################################################################################################################################### PAIR END MODE ###########################
if [ $PESE -eq 1 ]

echo Running in PE mode..........

then
        for i in $(cat $1);
                do
                array=(${i//,/ })
                echo -e 'Primer removal processing:\t'${array[0]}'\t'${array[1]}
                echo -e "${array[0]}" >> Stat.txt
                zcat ${array[0]} | echo $((`wc -l`/4)) >> Stat.txt
                echo -e "${array[1]}" >> Stat.txt
                zcat ${array[1]} | echo $((`wc -l`/4)) >> Stat.txt
                #Reverse Complement of the primers ##########################
                echo -e '>seq1\n'${array[3]} > 'primer.test'
                RcFwd=`revseq -sequence primer.test -stdout -auto | sed '1d'`
                echo -e '>seq1\n'${array[4]} > 'primer.test'
                RcRev=`revseq -sequence primer.test -stdout -auto | sed '1d'`
                rm primer.test
                #############################################################
                ## Fwd:${array[3]}, Rev:${array[4]}, Revcomp_Fwd:RcFwd, Revcomp_Rev:RcRev
                #TRIMMING OF TEH SEQUENCES BY LENGHT#####################
                if [ $TRIMBYCUT -eq 1 ]
                then
                        cutadapt -l $CUT_F ${array[0]} -o ${array[2]}.F.trimmed$CUT_F.fastq.gz
                        cutadapt -l $CUT_R ${array[1]} -o ${array[2]}.R.trimmed$CUT_R.fastq.gz
                        cutadapt -a ${array[3]}...$RcRev -A ${array[4]}...$RcFwd --discard-untrimmed -o  ${array[2]}_F_noadapt.fastq -p ${array[2]}_R_noadapt.fastq ${array[2]}.F.trimmed$CUT_F.fastq.gz ${array[2]}.R.trimmed$CUT_R.fastq.gz
                        cutadapt -a ${array[3]} -A ${array[4]} -o ${array[2]}_F_noadapt_2ndfilter.fastq  -p ${array[2]}_R_noadapt_2ndfilter.fastq ${array[2]}_F_noadapt.fastq ${array[2]}_R_noadapt.fastq
                        echo -e "PrimerRem\t${array[2]}_F_noadapt_2ndfilter.fastq\t$(count_fastq ${array[2]}_F_noadapt_2ndfilter.fastq)" >> Stat.txt
                        echo -e "PrimerRem\t${array[2]}_R_noadapt_2ndfilter.fastq\t$(count_fastq ${array[2]}_R_noadapt_2ndfilter.fastq)" >> Stat.txt
                        /home/dich/software/usearch_10.0.240_lin32 -fastq_mergepairs  ${array[2]}_F_noadapt_2ndfilter.fastq  -reverse ${array[2]}_R_noadapt_2ndfilter.fastq -fastqout ${array[2]}.merged.fq  -report ${array[2]}.merge_report -fastq_maxdiffs 0 -fastq_pctid $MERGEID -fastq_minqual 35 -relabel @ -tabbedout ${array[2]}.merge.map
                        /home/dich/software/usearch_10.0.240_lin32 -fastx_uniques ${array[2]}.merged.fq -fastaout ${array[2]}.uniq.fasta -sizeout -relabel ${array[2]}_
                        /home/dich/software/usearch_10.0.240_lin32 -fastq_mergepairs *_F_noadapt_2ndfilter.fastq -reverse *_R_noadapt_2ndfilter.fastq -fastqout all_merged.fq -relabel @ -report merge_usearch_rep.txt -fastq_maxdiffs 0 -fastq_pctid $MERGEID -fastq_minqual 35 -tabbedout merge_usearch_rep_tab.txt
                fi
                #TRIMMING BY TRIMMOMATIC #################################
                if [ $TRIMBYCUT -eq 0 ]
                then
                         cutadapt -a ${array[3]}...$RcRev -A ${array[4]}...$RcFwd --discard-untrimmed -o ${array[2]}_F_noadapt.fastq -p ${array[2]}_R_noadapt.fastq ${array[0]} ${array[1]}
                         java -jar /home/dich/software/Trimmomatic-0.38/trimmomatic-0.38.jar PE ${array[2]}_F_noadapt.fastq ${array[2]}_R_noadapt.fastq ${array[2]}_F_noadapt_Qfil_P.fq  ${array[2]}_F_noadapt_Qfil_UP.fq ${array[2]}_R_noadapt_Qfil_P.fq ${array[2]}_R_noadapt_Qfil_UP.fq ILLUMINACLIP:/home/dich/software/Trimmomatic-0.38/adapters/TruSeq2-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:$WINQ MINLEN:$MINL
                         echo -e "TrimmedbyTrimmomatic\t${array[2]}_F_noadapt_Qfil_P.fq\t$(count_fastq ${array[2]}_F_noadapt_Qfil_P.fq)" >> Stat.txt
                         echo -e "TrimmedbyTrimmomatic\t${array[2]}_R_noadapt_Qfil_P.fq\t$(count_fastq ${array[2]}_R_noadapt_Qfil_P.fq)" >> Stat.txt
                         /home/dich/software/usearch_10.0.240_lin32 -fastq_mergepairs ${array[2]}_F_noadapt_Qfil_P.fq -reverse ${array[2]}_R_noadapt_Qfil_P.fq -fastqout ${array[2]}.merged.fq -report ${array[2]}.merge_report -fastq_maxdiffs 0 -fastq_pctid $MERGEID -fastq_minqual 35 -relabel @ -tabbedout ${array[2]}.merge.map
                         /home/dich/software/usearch_10.0.240_lin32 -fastx_uniques ${array[2]}.merged.fq -fastaout ${array[2]}.uniq.fasta -sizeout -relabel ${array[2]}_
                        /home/dich/software/usearch_10.0.240_lin32 -fastq_mergepairs *F_noadapt_Qfil_P.fq -reverse *R_noadapt_Qfil_P.fq -fastqout all_merged.fq -relabel @ -report merge_usearch_rep.txt -fastq_maxdiffs 0 -fastq_pctid $MERGEID -fastq_minqual 35 -tabbedout merge_usearch_rep_tab.txt
                fi
                #/home/dich/scripts/needle2mut.pl ampli_ref.fa otu_less_00005.faa
                echo -e "Merge\t${array[2]}.merged.fq\t$(count_fastq ${array[2]}.merged.fq)" >> Stat.txt
                echo -e "Uniq\t${array[2]}.uniq.fasta\t$(count_fasta ${array[2]}.uniq.fasta)" >> Stat.txt
        done

        echo -e "#====Statitics all joint=====\n" >> Stat.txt
        echo filtering by file is finished... now will start the merge of the files, sorry that must be time consuming :/ 1. Vsearch 2. Usearch

#=====>> VSEARCH PIPE <=====================
cat *.merged.fq > all_cat_merged.fq
echo -e "#Vsearch##\nMerged_via_cat\tall_cat_merged.fq\t$(count_fastq all_cat_merged.fq)" >> Stat.txt
fastq_to_fasta -i all_cat_merged.fq -o all_cat_merged.fa
/home/dich/software/vsearch-2.8.1/bin/vsearch --derep_fulllength all_cat_merged.fa -sizeout -relabel var -output unique_from_cat_seqs.fa -threads 3
echo -e "Unique\tunique_from_cat_seqs.fa\t$(count_fasta unique_from_cat_seqs.fa)" >> Stat.txt
/home/dich/software/vsearch-2.8.1/bin/vsearch -usearch_global all_cat_merged.fa --db unique_from_cat_seqs.fa --id 1 --otutabout ASV_vsearch_counts.txt -threads 3
/home/dich/software/usearch_10.0.240_lin32 -otutab_trim ASV_vsearch_counts.txt -min_otu_freq 0.00005 -output ASV_vsearch_counts_00005_trim.txt
awk '{print $1}' ASV_vsearch_counts_00005_trim.txt | sed '1d' > id_ASV_vsearch_counts_00005_trim.list
/home/dich/software/usearch_10.0.240_lin32 -fastx_getseqs all_cat_merged.fq -labels id_ASV_vsearch_counts_00005_trim.list -fastaout id_ASV_vsearch_counts_00005_trim.fa
###########################################
echo done with Vsearch, we continue with Usearch
echo -e "#Usearch##\n" >> Stat.txt



#=====>> Usearch Arm  <=======================
#/home/dich/software/usearch_10.0.240_lin32 -fastq_mergepairs *_F_noadapt_2ndfilter.fastq -reverse *_R_noadapt_2ndfilter.fastq -fastqout all_merged.fq -relabel @ -report merge_usearch_rep.txt -fastq_maxdiffs 0 -fastq_pctid $MERGEID -fastq_minqual 35 -tabbedout merge_usearch_rep_tab.txt
echo -e "MergeDeNovo\tall_merged.fq\t$(count_fastq all_merged.fq)" >> Stat.txt
## To check the alignment options:  -tabbedout rep_tab  -alnout aln_out <- in the fastq_mergepairs
/home/dich/software/usearch_10.0.240_lin32 -fastx_uniques all_merged.fq -fastaout all.uniq.fasta -sizeout -relabel var
echo -e "Uniq\t$(count_fasta all.uniq.fasta)" >> Stat.txt
/home/dich/software/usearch_10.0.240_lin32 -otutab all_merged.fq -otus all.uniq.fasta -id 1 -otutabout otutab.txt -mapout map.txt -notmatched unmapped.fa -dbmatched otus_with_sizes.fa -sizeout
/home/dich/software/usearch_10.0.240_lin32 -otutab_trim otutab.txt -min_otu_freq 0.00005 -output otu_less_00005_trimmed.txt
/home/dich/software/usearch_10.0.240_lin32 -otutab_norm otu_less_00005_trimmed.txt -sample_size 300000 -output otu_less_00005_trimmed_norm.txt
awk '{print $1}' otu_less_00005_trimmed.txt | sed '1d' > lables_otu_less_00005.list
/home/dich/software/usearch_10.0.240_lin32 -fastx_getseqs all_merged.fq -labels lables_otu_less_00005.list -fastaout otu_less_00005.fa
######################################################
fi

echo Cool, everything is done.

if [ $CLEAN -eq 1 ]
        then
        rm *.trimmed$CUT_F.fastq.gz *.trimmed$CUT_R.fastq.gz *_2ndfilter.fastq *noadapt.fastq *merge.map *UP.fq *P.fq
        echo Clean..... Done.
fi

#============VARIANT description <<==================
#/home/dich/scripts/needle2mut.pl ampli_ref.fa otu_less_00005.faa
#####################################################
#rm *.fastq *.fa *.fq *.fasta *.faa *.txt *.json *.trimmed250.fastq.gz
