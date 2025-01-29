#$ -S /bin/bash
#$ -e /media/NewPsy/sge_logs/
#$ -o /media/NewPsy/sge_logs/

export FSLDIR=/usr/local/fsl
. /usr/local/fsl/etc/fslconf/fsl.sh

export PATH=$PATH:/opt/Shared_Software/R2016/bin
#export SPMDIR="/media/NewPsy/spm_mcr/spm_exec_v84"
export PATH="$PATH:/usr/lib64"
export SPM_LST=/opt/Shared_Software/R2016/toolbox/spm12/toolbox/LST
export ANTSPATH=/opt/Shared_Software/ants/install/bin
export PATH=${ANTSPATH}:${PATH}
export OMP_NUM_THREADS=1
export PATH=${FSLDIR}/bin:${PATH}
. $FSLDIR/etc/fslconf/fsl.sh


export FREESURFER_HOME=/opt/Shared_Software/FreeSurfer-7.3.2
source $FREESURFER_HOME/SetUpFreeSurfer.sh


# variabili di input
USER_NAME_ID=${1}
FileName_t13d=${2} # (zipped folder containig dicoms) OR (nii.gz)
#FileName_flair=${3}
FileName_dti=${3}
FileName_dti_corr=${4}
farmnp_t0=${5}
scl90_tot_t0=${6}
user_email=${7}

fileBaseName="${FileName_t13d%%.*}"
DRAMA_BPD_ID=$fileBaseName

# Usa la funzione `sed` per rimuovere i trattini e gli underscore
DRAMA_BPD_ID=$(echo "$DRAMA_BPD_ID" | sed 's/[-_]//g')


NP_PATH=/media/NewPsy

# Create Unique ID for Job
# Merge toghether date, time and 2 random numbers (YYYYMMDDhhmmssrr)
date_string=$(date '+%Y-%m-%d %H:%M:%S')
day=$(echo "${date_string}" | cut -d' ' -f1 | sed 's/-//g')
hour=$(echo "${date_string}" |cut -d' ' -f2 | sed 's/://g')
rnd=$(echo $RANDOM)
rnd=${rnd:0:2}

cd /media/NewPsy/sandbox


if [ ! -e ${USER_NAME_ID} ]; then
        mkdir ${USER_NAME_ID}
fi
cd ${USER_NAME_ID}
job_folder="${day}${hour}${rnd}"
mkdir $job_folder
cd $job_folder
job_home=$(pwd)

#mkdir /media/NewPsy/sandbox/${USER_NAME_ID}/${job_folder}
#cd /${job_folder}
mkdir /media/NewPsy/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}



#SET-UP VARIABLES
WD=/media/NewPsy/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}
mkdir ${WD}/T13D
mkdir ${WD}/DTI
mkdir ${WD}/FLAIR
mkdir ${WD}/RESULTS_DIR
mkdir $WD/PROCESSED_SUBJ
mkdir $WD/PROCESSED_SUBJ/FREESURFER
RESULTS_DIR=$WD/RESULTS_DIR
chmod 777 -R /media/NewPsy/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}

export SUBJECTS_DIR=$WD/PROCESSED_SUBJ/FREESURFER





                                                    ######################################## T13D ##############################################

if [ "${FileName_t13d}" != "NULL" ]; then
echo "T13D inserita"
#Controlla il tipo di file caricato T13D
if [ "${FileName_t13d##*.}" == "nii" ]; then
	mv /media/NewPsy/uploads/${FileName_t13d} ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/T13D/T13D_${DRAMA_BPD_ID}.nii
	gzip ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/T13D/T13D_${DRAMA_BPD_ID}.nii
elif [ "${FileName_t13d##*.}" == "gz" ]; then
	cp /media/NewPsy/uploads/${FileName_t13d} ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/T13D/T13D_${DRAMA_BPD_ID}.nii.gz #modified from mv to cp on 21/05/24
elif [ "${FileName_t13d##*.}" == "zip" ]; then
	echo "T13D in ingresso come DICOM"
	mkdir ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/T13D_o
	unzip -j /media/NewPsy/uploads/${FileName_t13d} -d ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/T13D_o

#Anonymization DICOMS
	if [ "${FileName_t13d##*.}" == "zip" ]; then
		#rm /media/NewPsy/uploads/${FileName_t13d}
		find ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/T13D_o -type f -follow -print > ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/list_dcm
		sed -i '/^$/d' ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/list_dcm

		chmod 777 -R ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}

		while read dcm_spect; do
			tipo=`dcmdump ${dcm_spect} | grep '(0002,0010)'`
			if [[ "$tipo" == *JPEG* ]]; then
				dcmdjpeg $dcm_spect  $dcm_spect
			fi
			dcmodify --no-backup -m "(0010,0010)= " ${dcm_spect}
			dcmodify --no-backup -m "(0010,0020)= " ${dcm_spect}
			dcmodify --no-backup -m "(0010,0030)= " ${dcm_spect}
			dcmodify --no-backup -m "(0010,0032)= " ${dcm_spect}
			dcmodify --no-backup -m "(0010,0050)= " ${dcm_spect}
			dcmodify --no-backup -m "(0010,1001)= " ${dcm_spect}
			dcmodify --no-backup -m "(0010,1040)= " ${dcm_spect}
			dcmodify --no-backup -m "(0010,1060)= " ${dcm_spect}
			dcmodify --no-backup -m "(0010,2154)= " ${dcm_spect}
			dcmodify --no-backup -m "(0038,0400)= " ${dcm_spect}
			dcmodify --no-backup -m "(0008,0020)= " ${dcm_spect}
			dcmodify --no-backup -m "(0008,0021)= " ${dcm_spect}
			dcmodify --no-backup -m "(0008,0022)= " ${dcm_spect}
			dcmodify --no-backup -m "(0008,0023)= " ${dcm_spect}
			dcmodify --no-backup -m "(0008,0030)= " ${dcm_spect}
			dcmodify --no-backup -m "(0008,0031)= " ${dcm_spect}
			dcmodify --no-backup -m "(0008,0032)= " ${dcm_spect}
			dcmodify --no-backup -m "(0008,0033)= " ${dcm_spect}
			dcmodify --no-backup -m "(0008,0050)= " ${dcm_spect}
			dcmodify --no-backup -m "(0008,0080)= " ${dcm_spect}
			dcmodify --no-backup -m "(0008,0081)= " ${dcm_spect}
			dcmodify --no-backup -m "(0008,0090)= " ${dcm_spect}
			dcmodify --no-backup -m "(0010,21b0)= " ${dcm_spect}
			dcmodify --no-backup -m "(0010,21c0)= " ${dcm_spect}
			dcmodify --no-backup -m "(0032,1032)= " ${dcm_spect}
			dcmodify --no-backup -m "(0032,1033)= " ${dcm_spect}
		done < ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/list_dcm
		rm ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/list_dcm
	fi

	cd ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/T13D_o
	mmv '* *' '#1#2'
	mmv '*(*' '#1#2'
	mmv '*)*' '#1#2'
	files=(*)
	tipo=`dcmdump ${files[0]} | grep '(0002,0010)'`
	if [[ "$tipo" == *JPEG* ]]; then

		ls > sub
		while read sub; do
			dcmdjpeg $sub  $sub
		done < sub
		rm sub
#cd $WD
	fi
	#dcm2niix
	export PATH=/opt/Shared_Software/dcm2niix/build/bin:$PATH
	dcm2niix ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/T13D_o
	img_t13d=$(ls ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/T13D_o/*.nii | head -1)
	mv ${img_t13d} ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/T13D/T13D_${DRAMA_BPD_ID}.nii
	gzip ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/T13D/T13D_${DRAMA_BPD_ID}.nii
fi


################## PRE-PROCESSING #####################

fslreorient2std ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/T13D/T13D_${DRAMA_BPD_ID}.nii.gz ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/T13D/T13D_${DRAMA_BPD_ID}_reor.nii.gz
rm ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/T13D/T13D_${DRAMA_BPD_ID}.nii.gz
mv ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/T13D/T13D_${DRAMA_BPD_ID}_reor.nii.gz ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/T13D/T13D_${DRAMA_BPD_ID}.nii.gz

#N4 correction
/opt/Shared_Software/ants/install/bin/N4BiasFieldCorrection -i ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/T13D/T13D_${DRAMA_BPD_ID}.nii.gz -o ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/T13D/T13D_${DRAMA_BPD_ID}_N4.nii.gz

#MNI registration
flirt -dof 12 -omat ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/T13D/T13D_${DRAMA_BPD_ID}_omat -in ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/T13D/T13D_${DRAMA_BPD_ID}_N4.nii.gz -ref /usr/local/fsl/data/standard/MNI152_T1_1mm.nii.gz -out ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/T13D/T13D_${DRAMA_BPD_ID}_N4_reg_to_MNI.nii.gz

####################### PROCESSING - FS ##################

recon-all -s ${DRAMA_BPD_ID} -i ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/T13D/T13D_${DRAMA_BPD_ID}_N4_reg_to_MNI.nii.gz -all
segmentThalamicNuclei.sh ${DRAMA_BPD_ID}
segmentHA_T1.sh ${DRAMA_BPD_ID} #segmentazione hippocampal subfields and nuclei of the amygdala
segmentBS.sh ${DRAMA_BPD_ID} #segmentazione Brainstem Substructures

# estrazione risultati

#recon-all
asegstats2table --subjects ${DRAMA_BPD_ID} --meas volume --all-segs -d comma --tablefile ${RESULTS_DIR}/${DRAMA_BPD_ID}_Cortical_Volume_aseg_stats.txt
aparcstats2table --subjects ${DRAMA_BPD_ID} --hemi lh --meas thickness -d comma --tablefile ${RESULTS_DIR}/${DRAMA_BPD_ID}_LH_thickness_aparc_stats.txt
aparcstats2table --subjects ${DRAMA_BPD_ID} --hemi rh --meas thickness -d comma --tablefile ${RESULTS_DIR}/${DRAMA_BPD_ID}_RH_thickness_aparc_stats.txt
paste -d, ${RESULTS_DIR}/${DRAMA_BPD_ID}_Cortical_Volume_aseg_stats.txt ${RESULTS_DIR}/${DRAMA_BPD_ID}_LH_thickness_aparc_stats.txt ${RESULTS_DIR}/${DRAMA_BPD_ID}_RH_thickness_aparc_stats.txt > ${RESULTS_DIR}/${DRAMA_BPD_ID}_Output_FS.csv 
#awk -F ',' -v OFS=',' '{print $117, $35, $25, $42, $19}' "${RESULTS_DIR}/${DRAMA_BPD_ID}_Output_FS.csv" > ${WD}/RESULTS_DIR/${DRAMA_BPD_ID}_Output_FS_features_selected.csv

#segmentHA_T1
sed 's/ /,/g' $WD/PROCESSED_SUBJ/FREESURFER/${DRAMA_BPD_ID}/mri/lh.hippoSfVolumes-T1.v22.txt > ${RESULTS_DIR}/lh_hippo.csv
csvtool transpose ${RESULTS_DIR}/lh_hippo.csv > ${RESULTS_DIR}/lh_hippo_transposed.csv
awk 'BEGIN{FS=OFS=","} NR==1{for(i=1;i<=NF;i++)$i="lh_"$i} 1' ${RESULTS_DIR}/lh_hippo_transposed.csv > ${RESULTS_DIR}/lh_hippo_transposed_prefixed.csv
sed 's/ /,/g' $WD/PROCESSED_SUBJ/FREESURFER/${DRAMA_BPD_ID}/mri/rh.hippoSfVolumes-T1.v22.txt > ${RESULTS_DIR}/rh_hippo.csv
csvtool transpose ${RESULTS_DIR}/rh_hippo.csv > ${RESULTS_DIR}/rh_hippo_transposed.csv
awk 'BEGIN{FS=OFS=","} NR==1{for(i=1;i<=NF;i++)$i="rh_"$i} 1' ${RESULTS_DIR}/rh_hippo_transposed.csv > ${RESULTS_DIR}/rh_hippo_transposed_prefixed.csv
sed 's/ /,/g' $WD/PROCESSED_SUBJ/FREESURFER/${DRAMA_BPD_ID}/mri/lh.amygNucVolumes-T1.v22.txt > ${RESULTS_DIR}/lh_amyg.csv
csvtool transpose ${RESULTS_DIR}/lh_amyg.csv > ${RESULTS_DIR}/lh_amyg_transposed.csv
awk 'BEGIN{FS=OFS=","} NR==1{for(i=1;i<=NF;i++)$i="lh_"$i} 1' ${RESULTS_DIR}/lh_amyg_transposed.csv > ${RESULTS_DIR}/lh_amyg_transposed_prefixed.csv
sed 's/ /,/g' $WD/PROCESSED_SUBJ/FREESURFER/${DRAMA_BPD_ID}/mri/rh.amygNucVolumes-T1.v22.txt > ${RESULTS_DIR}/rh_amyg.csv
csvtool transpose ${RESULTS_DIR}/rh_amyg.csv > ${RESULTS_DIR}/rh_amyg_transposed.csv
awk 'BEGIN{FS=OFS=","} NR==1{for(i=1;i<=NF;i++)$i="rh_"$i} 1' ${RESULTS_DIR}/rh_amyg_transposed.csv > ${RESULTS_DIR}/rh_amyg_transposed_prefixed.csv

#segmentBS
sed 's/ /,/g' $WD/PROCESSED_SUBJ/FREESURFER/${DRAMA_BPD_ID}/mri/brainstemSsVolumes.v13.txt > ${RESULTS_DIR}/Brain_Stem.csv
csvtool transpose ${RESULTS_DIR}/Brain_Stem.csv > ${RESULTS_DIR}/Brain_Stem_transposed.csv

#segmentThalamicNuclei
sed 's/ /,/g' $WD/PROCESSED_SUBJ/FREESURFER/${DRAMA_BPD_ID}/mri/ThalamicNuclei.v13.T1.volumes.txt > ${RESULTS_DIR}/thal_nuclei.csv
csvtool transpose ${RESULTS_DIR}/thal_nuclei.csv > ${RESULTS_DIR}/thal_nuclei_transposed.csv


paste -d, ${RESULTS_DIR}/${DRAMA_BPD_ID}_Output_FS.csv ${RESULTS_DIR}/lh_hippo_transposed_prefixed.csv ${RESULTS_DIR}/rh_hippo_transposed_prefixed.csv ${RESULTS_DIR}/lh_amyg_transposed_prefixed.csv ${RESULTS_DIR}/rh_amyg_transposed_prefixed.csv ${RESULTS_DIR}/Brain_Stem_transposed.csv ${RESULTS_DIR}/thal_nuclei_transposed.csv > ${RESULTS_DIR}/${DRAMA_BPD_ID}_Output_FS_TOTALE.csv

elif [ "${FileName_t13d}" == "NULL" ]; then

echo "T13D non presente"

fi

: <<'COMMENT'

                                ##########################################################FLAIR############################################################

if [ "${FileName_flair}" != "NULL" ]; then
echo "FLAIR inserita"

if [ "${FileName_flair##*.}" == "nii" ]; then
  	cp ${NP_PATH}/uploads/${FileName_flair} ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/FLAIR/FLAIR_${DRAMA_BPD_ID}.nii
  	gzip ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/FLAIR/FLAIR_${DRAMA_BPD_ID}.nii
elif [ "${FileName_flair##*.}" == "gz" ]; then
  	cp ${NP_PATH}/uploads/${FileName_flair} ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/FLAIR/FLAIR_${DRAMA_BPD_ID}.nii.gz
elif [ "${FileName_flair##*.}" == "zip" ]; then
  	mkdir ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/FLAIR_o
  	unzip -j ${NP_PATH}/uploads/${FileName_flair} -d ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/FLAIR_o
  	if [ "${FileName_flair##*.}" == "zip" ]; then
    		##anonymization DICOMS
    		rm ${NP_PATH}/uploads/${FileName_flair}
    		find ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/FLAIR_o -type f -follow -print > ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/list_dcm
    		ed -i '/^$/d' ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/list_dcm
    		while read dcm_spect; do
      			tipo=`dcmdump ${dcm_spect} | grep '(0002,0010)'`
      			if [[ "$tipo" == *JPEG* ]]; then
        			dcmdjpeg $dcm_spect $dcm_spect
      			fi
      			dcmodify --no-backup -m "(0010,0010)= " ${dcm_spect}
      			dcmodify --no-backup -m "(0010,0020)= " ${dcm_spect}
      			dcmodify --no-backup -m "(0010,0030)= " ${dcm_spect}
      			dcmodify --no-backup -m "(0010,0032)= " ${dcm_spect}
      			dcmodify --no-backup -m "(0010,0050)= " ${dcm_spect}
      			dcmodify --no-backup -m "(0010,1001)= " ${dcm_spect}
      			dcmodify --no-backup -m "(0010,1040)= " ${dcm_spect}
      			dcmodify --no-backup -m "(0010,1060)= " ${dcm_spect}
      			dcmodify --no-backup -m "(0010,2154)= " ${dcm_spect}
      			dcmodify --no-backup -m "(0038,0400)= " ${dcm_spect}
      			dcmodify --no-backup -m "(0008,0020)= " ${dcm_spect}
      			dcmodify --no-backup -m "(0008,0021)= " ${dcm_spect}
      			dcmodify --no-backup -m "(0008,0022)= " ${dcm_spect}
      			dcmodify --no-backup -m "(0008,0023)= " ${dcm_spect}
      			dcmodify --no-backup -m "(0008,0030)= " ${dcm_spect}
      			dcmodify --no-backup -m "(0008,0031)= " ${dcm_spect}
      			dcmodify --no-backup -m "(0008,0032)= " ${dcm_spect}
      			dcmodify --no-backup -m "(0008,0033)= " ${dcm_spect}
      			dcmodify --no-backup -m "(0008,0050)= " ${dcm_spect}
      			dcmodify --no-backup -m "(0008,0080)= " ${dcm_spect}
      			dcmodify --no-backup -m "(0008,0081)= " ${dcm_spect}
      			dcmodify --no-backup -m "(0008,0090)= " ${dcm_spect}
      			dcmodify --no-backup -m "(0010,21b0)= " ${dcm_spect}
      			dcmodify --no-backup -m "(0010,21c0)= " ${dcm_spect}
      			dcmodify --no-backup -m "(0032,1032)= " ${dcm_spect}
      			dcmodify --no-backup -m "(0032,1033)= " ${dcm_spect}
    			done < ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/list_dcm
    			rm ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/list_dcm
  	fi
  	cd ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/FLAIR_o
  	mmv '* *' '#1#2'
  	mmv '*(*' '#1#2'
  	mmv '*)*' '#1#2'
  	files=(*)
  	tipo=`dcmdump ${files[0]} | grep '(0002,0010)'`
  	if [[ "$tipo" == *JPEG* ]]; then
    		ls > sub
    		while read sub; do
     	 		dcmdjpeg $sub $sub
    		done < sub
    		rm sub
  	fi
  	export PATH=/opt/Shared_Software/dcm2niix/build/bin:$PATH
	dcm2niix ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/FLAIR_o
	img_flair=$(ls ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/FLAIR_o/*.nii | head -1)
	mv ${img_flair} ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/FLAIR/FLAIR_${DRAMA_BPD_ID}.nii
	gzip ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/FLAIR/FLAIR_${DRAMA_BPD_ID}.nii
fi

################## PRE-PROCESSING #####################

fslreorient2std ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/FLAIR/FLAIR_${DRAMA_BPD_ID}.nii.gz ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/FLAIR/FLAIR_${DRAMA_BPD_ID}_reor.nii.gz
rm ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/FLAIR/FLAIR_${DRAMA_BPD_ID}.nii.gz
mv ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/FLAIR/FLAIR_${DRAMA_BPD_ID}_reor.nii.gz ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/FLAIR/FLAIR_${DRAMA_BPD_ID}.nii.gz


#N4 correction
/opt/Shared_Software/ants/install/bin/N4BiasFieldCorrection -i ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/FLAIR/FLAIR_${DRAMA_BPD_ID}.nii.gz -o ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/FLAIR/FLAIR_${DRAMA_BPD_ID}_N4.nii.gz

#MNI registration
flirt -dof 12 -in ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/FLAIR/FLAIR_${DRAMA_BPD_ID}_N4.nii.gz -ref /usr/local/fsl/data/standard/MNI152_T1_1mm.nii.gz -out ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/FLAIR/FLAIR_${DRAMA_BPD_ID}_N4_reg_to_MNI.nii.gz

flirt -dof 12 -in ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/FLAIR/FLAIR_${DRAMA_BPD_ID}_N4_reg_to_MNI.nii.gz -ref ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/T13D/T13D_${DRAMA_BPD_ID}_N4_reg_to_MNI.nii.gz -out ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/FLAIR/FLAIR_${DRAMA_BPD_ID}_N4_reg_to_MNI_to_T1.nii.gz


##################### PROCESSING - LPA #######################

mkdir ${WD}/PROCESSED_SUBJ/LPA
mkdir ${WD}/PROCESSED_SUBJ/LPA/${DRAMA_BPD_ID}
chmod 777 ${WD}/PROCESSED_SUBJ/LPA
chmod 777 ${WD}/PROCESSED_SUBJ/LPA/${DRAMA_BPD_ID}

gunzip ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/FLAIR/FLAIR_${DRAMA_BPD_ID}_N4_reg_to_MNI_to_T1.nii.gz
gunzip ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/T13D/T13D_${DRAMA_BPD_ID}_N4_reg_to_MNI.nii.gz
cp ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/FLAIR/FLAIR_${DRAMA_BPD_ID}_N4_reg_to_MNI_to_T1.nii ${WD}/PROCESSED_SUBJ/LPA/${DRAMA_BPD_ID}
cp ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/T13D/T13D_${DRAMA_BPD_ID}_N4_reg_to_MNI.nii ${WD}/PROCESSED_SUBJ/LPA/${DRAMA_BPD_ID}

cd ${WD}/PROCESSED_SUBJ/LPA/${DRAMA_BPD_ID}
echo "Running LPA"

echo "
addpath('/opt/Shared_Software/R2016b/toolbox/spm12');
% List of open inputs
opengl software;
opengl('save','software');
nrun = 1; % enter the number of runs here
jobfile = {'${WD}/PROCESSED_SUBJ/LPA/${DRAMA_BPD_ID}/lpa_job.m'}; % Path to lpa_job.m
jobs = repmat(jobfile, 1, nrun);
inputs = cell(0, nrun);
for crun = 1:nrun
end
spm('defaults', 'PET');
spm_jobman('run', jobs, inputs{:});
quit()
" > ${WD}/PROCESSED_SUBJ/LPA/${DRAMA_BPD_ID}/lpa.m
  echo "
matlabbatch{1}.spm.tools.LST.lpa.data_F2 = {'${WD}/PROCESSED_SUBJ/LPA/${DRAMA_BPD_ID}/FLAIR_${DRAMA_BPD_ID}_N4_reg_to_MNI_to_T1.nii,1'};   %FLAIR
matlabbatch{1}.spm.tools.LST.lpa.data_coreg = {'${WD}/PROCESSED_SUBJ/LPA/${DRAMA_BPD_ID}/T13D_${DRAMA_BPD_ID}_N4_reg_to_MNI.nii,1'};   %T13D
matlabbatch{1}.spm.tools.LST.lpa.html_report = 1;
" > ${WD}/PROCESSED_SUBJ/LPA/${DRAMA_BPD_ID}/lpa_job.m
  /opt/Shared_Software/R2016b/bin/matlab -nodesktop -nosplash -r "lpa" -softwareopengl
  #cd ${WD}/PROCESSED_SUBJ/LPA/${DRAMA_BPD_ID}/
  lpa_html=$(find . -name "report_LST_lpa_*.html")
  vol=$(cat ${lpa_html}|sed '88!d' | sed -e 's/.*>\(.*\) ml.*/\1/')
  lesion=$(cat ${lpa_html}|sed '92!d' | sed -e 's/.*>\(.*\)<.*/\1/')
  echo "Lesion_volume","Lesion_number" >> ${RESULTS_DIR}/${DRAMA_BPD_ID}_Output_LPA.csv
  echo ${vol},${lesion} >> ${RESULTS_DIR}/${DRAMA_BPD_ID}_Output_LPA.csv

elif [ "${FileName_flair}" == "NULL" ]; then
echo "FLAIR non presente"
fi

COMMENT
                                                     ######################################### DTI #############################################
if [ "${FileName_dti}" != "NULL" ]; then
echo " DTI inserita"

if [ "${FileName_dti##*.}" == "zip" ]; then
    	mkdir ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI_o
    	mkdir ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI_o_corr
    	
    	                                                                            # DTI per correzione #
    if [ "${FileName_dti_corr}" != "NULL" ]; then
    		CORR=1
		echo "Presente DTI per correzione topup"
    		unzip -j ${NP_PATH}/uploads/${FileName_dti_corr} -d ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI_o_corr
    		
    		if [ $(find "${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI_o_corr" -type f -name "*.nii.gz" | wc -l) -eq 0 ]; then
        		#anonymization DICOMS
        		#rm ${NP_PATH}/Desktop/uploads/${FileName_dti}
        		find ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI_o_corr -type f -follow -print > ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/list_dcm
        		sed -i '/^$/d' ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/list_dcm
        		while read dcm_spect; do
            		tipo=`dcmdump ${dcm_spect} | grep '(0002,0010)'`
            		if [[ "$tipo" == *JPEG* ]]; then
                		dcmdjpeg $dcm_spect $dcm_spect
            		fi
            		dcmodify --no-backup -m "(0010,0010)= " ${dcm_spect}
            		dcmodify --no-backup -m "(0010,0020)= " ${dcm_spect}
            		dcmodify --no-backup -m "(0010,0030)= " ${dcm_spect}
            		dcmodify --no-backup -m "(0010,0032)= " ${dcm_spect}
            		dcmodify --no-backup -m "(0010,0050)= " ${dcm_spect}
            		dcmodify --no-backup -m "(0010,1001)= " ${dcm_spect}
            		dcmodify --no-backup -m "(0010,1040)= " ${dcm_spect}
            		dcmodify --no-backup -m "(0010,1060)= " ${dcm_spect}
            		dcmodify --no-backup -m "(0010,2154)= " ${dcm_spect}
            		dcmodify --no-backup -m "(0038,0400)= " ${dcm_spect}
            		dcmodify --no-backup -m "(0008,0020)= " ${dcm_spect}
            		dcmodify --no-backup -m "(0008,0021)= " ${dcm_spect}
            		dcmodify --no-backup -m "(0008,0022)= " ${dcm_spect}
            		dcmodify --no-backup -m "(0008,0023)= " ${dcm_spect}
            		dcmodify --no-backup -m "(0008,0030)= " ${dcm_spect}
            		dcmodify --no-backup -m "(0008,0031)= " ${dcm_spect}
            		dcmodify --no-backup -m "(0008,0032)= " ${dcm_spect}
            		dcmodify --no-backup -m "(0008,0033)= " ${dcm_spect}
            		dcmodify --no-backup -m "(0008,0050)= " ${dcm_spect}
            		dcmodify --no-backup -m "(0008,0080)= " ${dcm_spect}
            		dcmodify --no-backup -m "(0008,0081)= " ${dcm_spect}
            		dcmodify --no-backup -m "(0008,0090)= " ${dcm_spect}
            		dcmodify --no-backup -m "(0010,21b0)= " ${dcm_spect}
            		dcmodify --no-backup -m "(0010,21c0)= " ${dcm_spect}
            		dcmodify --no-backup -m "(0032,1032)= " ${dcm_spect}
            		dcmodify --no-backup -m "(0032,1033)= " ${dcm_spect}
        		done < ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/list_dcm
        		#rm ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/list_dcm
    		else
        		echo "\nNessun file DICOM trovato nella directory DTI_o per la DTI CORR\n."
    		fi
    		
		cd ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI_o_corr
    		ifexnii_corr1=$(find . -maxdepth 3 -type f -name "*.nii.gz" | head -1)
    		gunzip $ifexnii_corr1
    		ifexnii_corr=$(find . -maxdepth 3 -type f -name "*.nii" | head -1)
    		
    		if [ -z "$ifexnii_corr" ]; then
	    		echo "\nIn ingresso DICOM per la DTI CORR\n"
	    		
	    		
	      		files=$(find . -type f -follow -print | head -1)
	      		#echo "file dicom:"
	      		#echo $files
	      		tipo=`dcmdump ${files} | grep '(0002,0010)'`
	      		if [[ "$tipo" == *JPEG* ]]; then
				ls > sub
				while read sub; do
		  			dcmdjpeg $sub $sub
				done < sub
				rm sub
	      		fi
	      		dcm2niix ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI_o_corr
	      		img_dti=$(ls ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI_o_corr/*.nii | head -1)
	      		gzip ${img_dti}
	      		#rm ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI_o/co*.nii.gz
	      		img_dti=$(ls ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI_o_corr/*.nii.gz | head -1)
	      		if [ -z "$img_dti" ]; then
				cd ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI_o_corr
				img_dti=$(ls *.nii.gz | head -1)
				echo "IMMAGINI DTI"
				echo $img_dti
				mv ${img_dti} ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_FORCORR.nii.gz      
	      		else
				mv ${img_dti} ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_FORCORR.nii.gz      
	      		fi
	      		
	      		#Controlla se c'è il file bval dopo la conversione dcm2niix
	      		ifexbval_corr=$(find . -maxdepth 3 -type f -name "*.bval" | head -1)
	      		
	      		if [ -z "$ifexbval_corr" ]; then
	      			echo "\nbvals per DTI CORR non presenti\n"
	      			BVALSCORR=0
	      		else
	      			BVALSCORR=1
	      		
	      			img_bval=$(ls *.bval | head -1)
	      			img_bvec=$(ls *.bvec | head -1)
	      			linebval=$(cat $img_bval | sed '/^\s*$/d' | wc -l)
		      		
	      			echo "line_bval"
	      			echo ${linebval}
	      			if [ $linebval == '1' ]; then
	      		
	      				echo "check1"
					cat ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI_o_corr/${img_bval} | datamash -W transpose > ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_CORR.bvals
					cat ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI_o_corr/${img_bvec} | datamash -W transpose > ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_CORR.bvecs
	      			else
					echo "check2"
					mv $img_bval ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_CORR.bvals
					mv $img_bvec ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_CORR.bvecs
	      			fi
	      		
	      			grep "\S" ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_CORR.bvals > ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_CORR_nb.bvals
	      			rm ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_CORR.bvals
	      			mv ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_CORR_nb.bvals ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_CORR.bvals
	      			grep "\S" ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_CORR.bvecs > ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_CORR_nb.bvecs
	      			rm ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_CORR.bvecs
	      			mv ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_CORR_nb.bvecs ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_CORR.bvecs
	      			sed -i '/^\s*$/d' ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_CORR.bvecs
	      			sed -i '/^\s*$/d' ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_CORR.bvals
	      		fi
	    	else
    		
	    		echo "\nfile Nifti DTI FOR CORRECTION\n:"
	    		echo $ifexnii_corr
	    		mv ${ifexnii_corr} ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_FORCORR.nii
	      		gzip ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_FORCORR.nii
	      		
	      		#Controlla se è presente .bvals per DTI for correction
	      		ifexbvalscorr=$(find . -maxdepth 3 -type f -name "*.bvals" | head -1)
	      		
	      		if [ -z "$ifexbvalscorr" ]; then
	      			BVALSCORR=0
	      			echo "\nFile .bval non presente per DTI CORR\n"
	    			
	    		else
	    			BVALSCORR=1
	    			echo "\nFile .bval presente per DTI CORR\n"
	    			img_bvalcorr=$(ls *.bval | head -1)
	    			linebvalcorr=$(cat $img_bvalcorr | sed '/^\s*$/d' | wc -l)
	    			
	    			if [ $linebvalcorr == '1' ]; then
	      		
	      				cat ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI_o_corr/${img_bvalcorr} | datamash -W transpose > ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_CORR.bvals
					#cat ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI_o/${img_bvec} | datamash -W transpose > ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvecs
	      			else
					
					mv $img_bvalcorr ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_CORR.bvals
					#mv $img_bvec ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvecs
	      			fi
	      			grep "\S" ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_CORR.bvals > ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_CORR_nb.bvals
	      			rm ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_CORR.bvals
	      			mv ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_CORR_nb.bvals ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_CORR.bvals
	      			#grep "\S" ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvecs > ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_nb.bvecs
	      			#rm ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvecs
	      			#mv ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_nb.bvecs ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvecs
	      			#sed -i '/^\s*$/d' ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvecs
	      			sed -i '/^\s*$/d' ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_CORR.bvals
	    		fi
		fi
		
    elif [ "${FileName_dti_corr}" == "NULL" ]; then	
		CORR=0
    fi

    	unzip -j ${NP_PATH}/uploads/${FileName_dti} -d ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI_o
    	
		                                                                  # DTI reale #
	
	if [ "${FileName_dti##*.}" == "zip" ]; then
		#rm ${NP_PATH}/uploads/${FileName_dti}
    		# Controlla se ci sono file .dcm nella directory DTI_o
    		if [ $(find "${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI_o" -type f -name "*.nii.gz" | wc -l) -eq 0 ]; then
        		#anonymization DICOMS
        		#rm ${NP_PATH}/uploads/${FileName_dti}
        		find ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI_o -type f -follow -print > ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/list_dcm
        		sed -i '/^$/d' ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/list_dcm
        		while read dcm_spect; do
            		tipo=`dcmdump ${dcm_spect} | grep '(0002,0010)'`
            		if [[ "$tipo" == *JPEG* ]]; then
                		dcmdjpeg $dcm_spect $dcm_spect
            		fi
            		dcmodify --no-backup -m "(0010,0010)= " ${dcm_spect}
            		dcmodify --no-backup -m "(0010,0020)= " ${dcm_spect}
            		dcmodify --no-backup -m "(0010,0030)= " ${dcm_spect}
            		dcmodify --no-backup -m "(0010,0032)= " ${dcm_spect}
            		dcmodify --no-backup -m "(0010,0050)= " ${dcm_spect}
            		dcmodify --no-backup -m "(0010,1001)= " ${dcm_spect}
            		dcmodify --no-backup -m "(0010,1040)= " ${dcm_spect}
            		dcmodify --no-backup -m "(0010,1060)= " ${dcm_spect}
            		dcmodify --no-backup -m "(0010,2154)= " ${dcm_spect}
            		dcmodify --no-backup -m "(0038,0400)= " ${dcm_spect}
            		dcmodify --no-backup -m "(0008,0020)= " ${dcm_spect}
            		dcmodify --no-backup -m "(0008,0021)= " ${dcm_spect}
            		dcmodify --no-backup -m "(0008,0022)= " ${dcm_spect}
            		dcmodify --no-backup -m "(0008,0023)= " ${dcm_spect}
            		dcmodify --no-backup -m "(0008,0030)= " ${dcm_spect}
            		dcmodify --no-backup -m "(0008,0031)= " ${dcm_spect}
            		dcmodify --no-backup -m "(0008,0032)= " ${dcm_spect}
            		dcmodify --no-backup -m "(0008,0033)= " ${dcm_spect}
            		dcmodify --no-backup -m "(0008,0050)= " ${dcm_spect}
            		dcmodify --no-backup -m "(0008,0080)= " ${dcm_spect}
            		dcmodify --no-backup -m "(0008,0081)= " ${dcm_spect}
            		dcmodify --no-backup -m "(0008,0090)= " ${dcm_spect}
            		dcmodify --no-backup -m "(0010,21b0)= " ${dcm_spect}
            		dcmodify --no-backup -m "(0010,21c0)= " ${dcm_spect}
            		dcmodify --no-backup -m "(0032,1032)= " ${dcm_spect}
            		dcmodify --no-backup -m "(0032,1033)= " ${dcm_spect}
        		done < ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/list_dcm
        		#rm ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/list_dcm
    		else
        		echo "\nNessun file DICOM trovato nella directory DTI_o\n"
    		fi
	fi

    	cd ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI_o
    	chmod 777 -R ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}
    	ifexnii1=$(find . -maxdepth 3 -type f -name "*.nii.gz" | head -1)
    	gunzip $ifexnii1
    	ifexnii=$(find . -maxdepth 3 -type f -name "*.nii" | head -1)
    	echo "\nfile Nifti:\n"
    	echo $ifexnii
    	cd ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI_o
    	mmv '* *' '#1#2'
   	mmv '*(*' '#1#2'
    	mmv '*)*' '#1#2'
    	echo "check"
    	if [ -z "$ifexnii" ]; then
    		echo "\nIn ingresso DICOM\n"
    		JSON=1
    		
      		files=$(find . -type f -follow -print | head -1)
      		echo "file dicom:"
      		echo $files
      		tipo=`dcmdump ${files} | grep '(0002,0010)'`
      		if [[ "$tipo" == *JPEG* ]]; then
        		ls > sub
        		while read sub; do
          			dcmdjpeg $sub $sub
        		done < sub
        		rm sub
      		fi
      		dcm2niix ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI_o
      		img_dti=$(ls ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI_o/*.nii | head -1)
      		gzip ${img_dti}
      		#rm ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI_o/co*.nii.gz
      		img_dti=$(ls ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI_o/*.nii.gz | head -1)
      		if [ -z "$img_dti" ]; then
        		cd ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI_o
        		img_dti=$(ls *.nii.gz | head -1)
        		echo "IMMAGINI DTI"
        		echo $img_dti
        		mv ${img_dti} ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.nii.gz      
      		else
        		mv ${img_dti} ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.nii.gz      
      		fi
      		img_bval=$(ls *.bval | head -1)
      		img_bvec=$(ls *.bvec | head -1)
      		linebval=$(cat $img_bval | sed '/^\s*$/d' | wc -l)
      		
      		#check for PED in .json
      		PE=$(jq -r '.PhaseEncodingDirection' *.json)
		if [[ $PE == "j-" ]]; then
    			PE="AP"
		elif [[ $PE == "j" ]]; then
    			PE="PA"
    		else
    			PE="AP"
		fi
		echo "\nPhaseEncodingDirection: ${PE}\n"
		
      		#echo "line_bval"
      		echo ${linebval}
      		if [ $linebval == '1' ]; then
      		
      			echo "check1"
        		cat ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI_o/${img_bval} | datamash -W transpose > ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvals
        		cat ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI_o/${img_bvec} | datamash -W transpose > ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvecs
      		else
        		echo "check2"
        		mv $img_bval ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvals
        		mv $img_bvec ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvecs
      		fi
      		grep "\S" ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvals > ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_nb.bvals
      		rm ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvals
      		mv ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_nb.bvals ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvals
      		grep "\S" ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvecs > ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_nb.bvecs
      		rm ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvecs
      		mv ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_nb.bvecs ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvecs
      		sed -i '/^\s*$/d' ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvecs
      		sed -i '/^\s*$/d' ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvals
    	else
      		echo "\nIn ingresso Nifti\n"
      		ifexjson=$(find . -maxdepth 3 -type f -name "*.json")
      		if [ -z "$ifexjson" ]; then
      			JSON=0
      			echo "\nFile .json non presente: PE_dir: AP\n"
    			
    		else
    			JSON=1
    			echo "\nRecupero PE_dir da file .json\n"
    			
    			PE=$(jq -r '.PhaseEncodingDirection' *.json)
			if [[ $PE == "j-" ]]; then
    				PE="AP"
			elif [[ $PE == "j" ]]; then
    				PE="PA"
    			else
    				PE="AP"
			fi
		fi
		echo "\nPhaseEncodingDirection: ${PE}\n"
		
      		mv ${ifexnii} ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.nii
      		gzip ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.nii
      		img_bval=$(ls *.bval | head -1)
      		img_bvec=$(ls *.bvec | head -1)
      		linebval=$(cat $img_bval | sed '/^\s*$/d' | wc -l)
      		#echo "line_bval"
      		echo ${linebval}
      		if [ $linebval == '1' ]; then
      		
      			echo "\nTraspongo bvals e bvecs\n"
        		cat ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI_o/${img_bval} | datamash -W transpose > ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvals
        		cat ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI_o/${img_bvec} | datamash -W transpose > ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvecs
      		else
        		echo "\nNON Traspongo bvals e bvecs\n"
        		mv $img_bval ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvals
        		mv $img_bvec ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvecs
      		fi
      		grep "\S" ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvals > ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_nb.bvals
      		rm ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvals
      		mv ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_nb.bvals ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvals
      		grep "\S" ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvecs > ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_nb.bvecs
      		rm ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvecs
      		mv ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_nb.bvecs ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvecs
      		sed -i '/^\s*$/d' ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvecs
      		sed -i '/^\s*$/d' ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvals
    	fi
    	cd $WD
fi


#################### PRE-PROCESSING ###################

echo "\n\n########## START OF DTI PREPROCESSING ##########\n\n"

fslreorient2std ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.nii.gz ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_reor.nii.gz
rm ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.nii.gz
mv ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_reor.nii.gz ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.nii.gz


#CORREZIONE EDDY PER LA DWI

if [ ${CORR} == '1' ]; then 
	if [ ${BVALSCORR} == '1' ]; then 

		#esegui correzione topup
		#echo $CORR
		echo "\n\nCORREZIONE TOP-UP e BVAL CORR presente\n\n"
		
		fslreorient2std ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_FORCORR.nii.gz ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_reor_corr.nii.gz
		rm ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_FORCORR.nii.gz
		mv ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_reor_corr.nii.gz ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_FORCORR.nii.gz

		if [ ${JSON} == '1' ]; then 
			echo "\nEseguo TOPUP con PE_dir presa dal .json\n"

			# DTI per correzione 

			/opt/Shared_Software/mrtrix3/bin/dwidenoise ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_FORCORR.nii.gz  ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_FORCORR_denoise.nii.gz
			/opt/Shared_Software/mrtrix3/bin/mrdegibbs ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_FORCORR_denoise.nii.gz ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_FORCORR_denoise_degibbs.nii.gz
			
					
			# Trova indici con bval = 0
			indices_corr=($(awk '$1 == 0 {print NR}' ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_CORR.bvals))
			i=1
			for index in "${indices_corr[@]}"
			do
	  			echo "Indice della riga con valore 0: $index"
	  			echo ${i}
	  			
	  			#estrai volume con bval = 0
				fslroi ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_FORCORR_denoise_degibbs.nii.gz ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_FORCORR_denoise_degibbs_1vol_${i}.nii.gz $((${index}-1)) 1 
				i=$((i + 1))
			done

			# Costruisci il comando fslmaths dinamicamente
			fslmaths_command="fslmaths"

			# Aggiungi ciascun file DTI al comando fslmaths
			for ((j=1; j<i; j++))
			do
	  			fslmaths_command+=" ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_FORCORR_denoise_degibbs_1vol_${j}.nii.gz -add"
			done

			# Rimuove l'ultima occorrenza di '-add'
			fslmaths_command="${fslmaths_command% -add}"

			# Calcola la divisione dinamica per il numero di file DTI
			fslmaths_command+=" -div $((${i}-1)) ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_FORCORR_denoise_degibbs_MEDIA_vols.nii.gz"

			# Esegui il comando fslmaths
			echo "Eseguendo comando: $fslmaths_command"
			eval $fslmaths_command


			# DTI reale
					
			/opt/Shared_Software/mrtrix3/bin/dwidenoise ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.nii.gz  ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise.nii.gz
			/opt/Shared_Software/mrtrix3/bin/mrdegibbs ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise.nii.gz ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs.nii.gz
			
			# Trova indici con bval = 0
			indices_real=($(awk '$1 == 0 {print NR}' ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvals))
			i=1
			for index in "${indices_real[@]}"
			do
	  			echo "Indice della riga con valore 0: $index"
	  			echo ${i}
	  			
	  			#estrai volume con bval = 0
				fslroi ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs.nii.gz ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs_1vol_${i}.nii.gz $((${index}-1)) 1 
				i=$((i + 1))
			done

			# Costruisci il comando fslmaths dinamicamente
			fslmaths_command="fslmaths"

			# Aggiungi ciascun file DTI al comando fslmaths
			for ((j=1; j<i; j++))
			do
	  			fslmaths_command+=" ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs_1vol_${j}.nii.gz -add"
			done

			# Rimuove l'ultima occorrenza di '-add'
			fslmaths_command="${fslmaths_command% -add}"

			# Calcola la divisione dinamica per il numero di file DTI
			fslmaths_command+=" -div $((${i}-1)) ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs_MEDIA_vols.nii.gz"

			# Esegui il comando fslmaths
			echo "Eseguendo comando: $fslmaths_command"
			eval $fslmaths_command
			
			
			#merge fra le due medie 
			fslmerge -t  ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_MEDIATED.nii.gz  ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs_MEDIA_vols.nii.gz ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_FORCORR_denoise_degibbs_MEDIA_vols.nii.gz 


			/opt/Shared_Software/mrtrix3/bin/dwifslpreproc  ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs.nii.gz ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs_PREPROC.nii.gz -rpe_pair -fslgrad  "${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvecs" "${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvals" -se_epi ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_MEDIATED.nii.gz -pe_dir ${PE} -eddy_options " --slm=linear"

			/opt/Shared_Software/mrtrix3/bin/dwibiascorrect ants ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs_PREPROC.nii.gz ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs_PREPROC_BIASC.nii.gz -fslgrad "${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvecs" "${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvals"



		elif [ ${JSON} == '0' ]; then
			echo "\nEseguo TOPUP con PE_dir imposta\n"

			# DTI per correzione 

			/opt/Shared_Software/mrtrix3/bin/dwidenoise ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_FORCORR.nii.gz  ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_FORCORR_denoise.nii.gz
			/opt/Shared_Software/mrtrix3/bin/mrdegibbs ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_FORCORR_denoise.nii.gz ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_FORCORR_denoise_degibbs.nii.gz

			# Trova indici con bval = 0
			indices_corr=($(awk '$1 == 0 {print NR}' ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_CORR.bvals))
			i=1
			for index in "${indices_corr[@]}"
			do
	  			echo "Indice della riga con valore 0: $index"
	  			echo ${i}
	  			
	  			#estrai volume con bval = 0
				fslroi ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_FORCORR_denoise_degibbs.nii.gz ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_FORCORR_denoise_degibbs_1vol_${i}.nii.gz $((${index}-1)) 1 
				i=$((i + 1))
			done

			# Costruisci il comando fslmaths dinamicamente
			fslmaths_command="fslmaths"

			# Aggiungi ciascun file DTI al comando fslmaths
			for ((j=1; j<i; j++))
			do
	  			fslmaths_command+=" ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_FORCORR_denoise_degibbs_1vol_${j}.nii.gz -add"
			done

			# Rimuove l'ultima occorrenza di '-add'
			fslmaths_command="${fslmaths_command% -add}"

			# Calcola la divisione dinamica per il numero di file DTI
			fslmaths_command+=" -div $((${i}-1)) ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_FORCORR_denoise_degibbs_MEDIA_vols.nii.gz"

			# Esegui il comando fslmaths
			echo "Eseguendo comando: $fslmaths_command"
			eval $fslmaths_command

			# DTI reale

			/opt/Shared_Software/mrtrix3/bin/dwidenoise ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.nii.gz  ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise.nii.gz
			/opt/Shared_Software/mrtrix3/bin/mrdegibbs ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise.nii.gz ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs.nii.gz

			# Trova indici con bval = 0
			indices_real=($(awk '$1 == 0 {print NR}' ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvals))
			i=1
			for index in "${indices_real[@]}"
			do
	  			echo "Indice della riga con valore 0: $index"
	  			echo ${i}
	  			
	  			#estrai volume con bval = 0
				fslroi ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs.nii.gz ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs_1vol_${i}.nii.gz $((${index}-1)) 1 
				i=$((i + 1))
			done

			# Costruisci il comando fslmaths dinamicamente
			fslmaths_command="fslmaths"

			# Aggiungi ciascun file DTI al comando fslmaths
			for ((j=1; j<i; j++))
			do
	  			fslmaths_command+=" ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs_1vol_${j}.nii.gz -add"
			done

			# Rimuove l'ultima occorrenza di '-add'
			fslmaths_command="${fslmaths_command% -add}"

			# Calcola la divisione dinamica per il numero di file DTI
			fslmaths_command+=" -div $((${i}-1)) ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs_MEDIA_vols.nii.gz"

			# Esegui il comando fslmaths
			echo "Eseguendo comando: $fslmaths_command"
			eval $fslmaths_command

			#merge fra le due medie 
			fslmerge -t  ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_MEDIATED.nii.gz  ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs_MEDIA_vols.nii.gz ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_FORCORR_denoise_degibbs_MEDIA_vols.nii.gz 


			/opt/Shared_Software/mrtrix3/bin/dwifslpreproc  ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs.nii.gz ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs_PREPROC.nii.gz -rpe_pair -fslgrad  "${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvecs" "${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvals" -se_epi ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_MEDIATED.nii.gz -pe_dir AP -eddy_options " --slm=linear"

			/opt/Shared_Software/mrtrix3/bin/dwibiascorrect ants ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs_PREPROC.nii.gz ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs_PREPROC_BIASC.nii.gz -fslgrad "${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvecs" "${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvals"

		fi  

	elif [ ${BVALSCORR} == '0' ]; then 

		echo "\n\nCORREZIONE TOP-UP e BVAL CORR NON presente\n\n"
		
		#esegui correzione topup
		fslreorient2std ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_FORCORR.nii.gz ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_reor_corr.nii.gz
		rm ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_FORCORR.nii.gz
		mv ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_reor_corr.nii.gz ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_FORCORR.nii.gz
		

		if [ ${JSON} == '1' ]; then 
			echo "\nEseguo TOPUP con PE_dir presa dal .json\n"

			# DTI per correzione 

			/opt/Shared_Software/mrtrix3/bin/dwidenoise ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_FORCORR.nii.gz  ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_FORCORR_denoise.nii.gz
			/opt/Shared_Software/mrtrix3/bin/mrdegibbs ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_FORCORR_denoise.nii.gz ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_FORCORR_denoise_degibbs.nii.gz
			
			#TROVA I VOLUMI CON BVAL = 0 SENZA AVERE I BVALS
			
			# Ottieni il numero di volumi con fslinfo
			dim4_value=$(fslinfo ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_FORCORR.nii.gz | grep ^dim4 | awk '{print $2}')
			echo "Il numero di volumi della DTI FOR CORR è: $dim4_value"

			# Estrai ogni volume dalla FORCORR originale usando fslroi
			for ((j=0; j<$dim4_value; j++))
			do
			    fslroi ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_FORCORR.nii.gz ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_FORCORR_vol_${j}_o.nii.gz ${j} 1
			done 
			
			cd ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/
			# Inizializza il contatore e il file di output
			
			i=0
			output_file="volume_intensities.txt"

			# Rimuovi file di output esistenti
			rm -f $output_file

			# Calcola l'intensità media per ogni volume
			while [ -e "DTI_${DRAMA_BPD_ID}_FORCORR_vol_${i}_o.nii.gz" ]; do
			  mean_intensity=$(mrstats "DTI_${DRAMA_BPD_ID}_FORCORR_vol_${i}_o.nii.gz" -output mean)
			  echo "Volume ${i} ${mean_intensity}" >> $output_file
			  i=$((i + 1))
			done

			# Visualizza le intensità
			cat $output_file

			# File contenente le intensità
			file="volume_intensities.txt"

			# Calcola la media delle intensità
			mean=$(awk '{sum += $3; count++} END {print sum/count}' $file)

			# Calcola la deviazione standard delle intensità
			std_dev=$(awk -v mean="$mean" '{sum += ($3 - mean)^2; count++} END {print sqrt(sum/count)}' $file)

			echo "Media delle intensità: $mean"
			echo "Deviazione standard: $std_dev"

			# Inizializza una variabile per salvare i volumi con intensità alta
			high_intensity_volumes="bval0_volumes_FORCORR.txt"
			rm -f $high_intensity_volumes

			# Variabile per tenere traccia se è stato trovato almeno un volume con intensità > mean + 2 * std_dev
			found_high_intensity=false

			# Legge il file riga per riga
			while IFS= read -r line; do
			    volume=$(echo "$line" | awk '{print $2}')
			    intensity=$(echo "$line" | awk '{print $3}')

			    # Verifica se l'intensità è maggiore di mean + 2 * std_dev
			    if (( $(echo "$intensity > $mean + 2 * $std_dev" | bc -l) )); then
				echo "$volume" >> $high_intensity_volumes
				found_high_intensity=true
			    fi
			done < $file

			# Se nessun volume ha intensità > mean + 2 * std_dev, aggiunge tutti i volumi al file
			if [ "$found_high_intensity" = false ]; then
			    while IFS= read -r line; do
				volume=$(echo "$line" | awk '{print $2}')
				echo "$volume" >> $high_intensity_volumes
			    done < $file
			fi
			
			cat $high_intensity_volumes

			file_bval0_vol="bval0_volumes_FORCORR.txt"

			# Leggi ogni riga del file e estrae il valore numerico
			i=1
			for value in $(cat "$file_bval0_vol"); do
			    	echo "Volume con bval = 0: $value"
			    	
			    	#estrai volume con bval = 0
				fslroi ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_FORCORR_denoise_degibbs.nii.gz ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_FORCORR_denoise_degibbs_1vol_${i}.nii.gz $value 1
				i=$((i + 1))
			done
			
			# Costruisci il comando fslmaths dinamicamente
			fslmaths_command="fslmaths"

			# Aggiungi ciascun file DTI al comando fslmaths
			for ((j=1; j<i; j++))
			do
	  			fslmaths_command+=" ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_FORCORR_denoise_degibbs_1vol_${j}.nii.gz -add"
			done

			# Rimuove l'ultima occorrenza di '-add'
			fslmaths_command="${fslmaths_command% -add}"

			# Calcola la divisione dinamica per il numero di file DTI
			fslmaths_command+=" -div $((${i}-1)) ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_FORCORR_denoise_degibbs_MEDIA_vols.nii.gz"

			# Esegui il comando fslmaths
			echo "Eseguendo comando: $fslmaths_command"
			eval $fslmaths_command		
		


			# DTI reale
					
			/opt/Shared_Software/mrtrix3/bin/dwidenoise ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.nii.gz  ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise.nii.gz
			/opt/Shared_Software/mrtrix3/bin/mrdegibbs ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise.nii.gz ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs.nii.gz
			
			# Trova indici con bval = 0
			indices_real=($(awk '$1 == 0 {print NR}' ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvals))
			i=1
			for index in "${indices_real[@]}"
			do
	  			echo "Indice della riga con valore 0: $index"
	  			echo ${i}
	  			
	  			#estrai volume con bval = 0
				fslroi ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs.nii.gz ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs_1vol_${i}.nii.gz $((${index}-1)) 1 
				i=$((i + 1))
			done

			# Costruisci il comando fslmaths dinamicamente
			fslmaths_command="fslmaths"

			# Aggiungi ciascun file DTI al comando fslmaths
			for ((j=1; j<i; j++))
			do
	  			fslmaths_command+=" ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs_1vol_${j}.nii.gz -add"
			done

			# Rimuove l'ultima occorrenza di '-add'
			fslmaths_command="${fslmaths_command% -add}"

			# Calcola la divisione dinamica per il numero di file DTI
			fslmaths_command+=" -div $((${i}-1)) ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs_MEDIA_vols.nii.gz"

			# Esegui il comando fslmaths
			echo "Eseguendo comando: $fslmaths_command"
			eval $fslmaths_command
			
			
			#merge fra le due medie 
			fslmerge -t  ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_MEDIATED.nii.gz  ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs_MEDIA_vols.nii.gz ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_FORCORR_denoise_degibbs_MEDIA_vols.nii.gz 


			/opt/Shared_Software/mrtrix3/bin/dwifslpreproc  ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs.nii.gz ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs_PREPROC.nii.gz -rpe_pair -fslgrad  "${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvecs" "${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvals" -se_epi ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_MEDIATED.nii.gz -pe_dir ${PE} -eddy_options " --slm=linear"

			/opt/Shared_Software/mrtrix3/bin/dwibiascorrect ants ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs_PREPROC.nii.gz ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs_PREPROC_BIASC.nii.gz -fslgrad "${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvecs" "${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvals"



		elif [ ${JSON} == '0' ]; then
			echo "\nEseguo TOPUP con PE_dir imposta\n"

			# DTI per correzione 

			/opt/Shared_Software/mrtrix3/bin/dwidenoise ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_FORCORR.nii.gz  ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_FORCORR_denoise.nii.gz
	/opt/Shared_Software/mrtrix3/bin/mrdegibbs ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_FORCORR_denoise.nii.gz ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_FORCORR_denoise_degibbs.nii.gz
	
			#TROVA I VOLUMI CON BVAL = 0 SENZA AVERE I BVALS
			
			# Ottieni il numero di volumi con fslinfo
			dim4_value=$(fslinfo ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_FORCORR.nii.gz | grep ^dim4 | awk '{print $2}')
			echo "Il numero di volumi della DTI FOR CORR è: $dim4_value"

			# Estrai ogni volume dalla FORCORR originale usando fslroi
			for ((j=0; j<$dim4_value; j++))
			do
			    fslroi ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_FORCORR.nii.gz ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_FORCORR_vol_${j}_o.nii.gz ${j} 1
			done 
			
			cd ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/
			# Inizializza il contatore e il file di output
			
			i=0
			output_file="volume_intensities.txt"

			# Rimuovi file di output esistenti
			rm -f $output_file

			# Calcola l'intensità media per ogni volume
			while [ -e "DTI_${DRAMA_BPD_ID}_FORCORR_vol_${i}_o.nii.gz" ]; do
			  mean_intensity=$(mrstats "DTI_${DRAMA_BPD_ID}_FORCORR_vol_${i}_o.nii.gz" -output mean)
			  echo "Volume ${i} ${mean_intensity}" >> $output_file
			  i=$((i + 1))
			done

			# Visualizza le intensità
			cat $output_file

			# File contenente le intensità
			file="volume_intensities.txt"

			# Calcola la media delle intensità
			mean=$(awk '{sum += $3; count++} END {print sum/count}' $file)

			# Calcola la deviazione standard delle intensità
			std_dev=$(awk -v mean="$mean" '{sum += ($3 - mean)^2; count++} END {print sqrt(sum/count)}' $file)

			echo "Media delle intensità: $mean"
			echo "Deviazione standard: $std_dev"

			# Inizializza una variabile per salvare i volumi con intensità alta
			high_intensity_volumes="bval0_volumes_FORCORR.txt"
			rm -f $high_intensity_volumes

			# Variabile per tenere traccia se è stato trovato almeno un volume con intensità > mean + 2 * std_dev
			found_high_intensity=false

			# Legge il file riga per riga
			while IFS= read -r line; do
			    volume=$(echo "$line" | awk '{print $2}')
			    intensity=$(echo "$line" | awk '{print $3}')

			    # Verifica se l'intensità è maggiore di mean + 2 * std_dev
			    if (( $(echo "$intensity > $mean + 2 * $std_dev" | bc -l) )); then
				echo "$volume" >> $high_intensity_volumes
				found_high_intensity=true
			    fi
			done < $file

			# Se nessun volume ha intensità > mean + 2 * std_dev, aggiunge tutti i volumi al file
			if [ "$found_high_intensity" = false ]; then
			    while IFS= read -r line; do
				volume=$(echo "$line" | awk '{print $2}')
				echo "$volume" >> $high_intensity_volumes
			    done < $file
			fi
			
			cat $high_intensity_volumes

			file_bval0_vol="bval0_volumes_FORCORR.txt"

			# Leggi ogni riga del file e estrae il valore numerico
			i=1
			for value in $(cat "$file_bval0_vol"); do
			    	echo "Volume con bval = 0: $value"
			    	
			    	#estrai volume con bval = 0
				fslroi ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_FORCORR_denoise_degibbs.nii.gz ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_FORCORR_denoise_degibbs_1vol_${i}.nii.gz $value 1
				i=$((i + 1))
			done
			
			# Costruisci il comando fslmaths dinamicamente
			fslmaths_command="fslmaths"

			# Aggiungi ciascun file DTI al comando fslmaths
			for ((j=1; j<i; j++))
			do
	  			fslmaths_command+=" ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_FORCORR_denoise_degibbs_1vol_${j}.nii.gz -add"
			done

			# Rimuove l'ultima occorrenza di '-add'
			fslmaths_command="${fslmaths_command% -add}"

			# Calcola la divisione dinamica per il numero di file DTI
			fslmaths_command+=" -div $((${i}-1)) ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_FORCORR_denoise_degibbs_MEDIA_vols.nii.gz"

			# Esegui il comando fslmaths
			echo "Eseguendo comando: $fslmaths_command"
			eval $fslmaths_command		
				
			# DTI reale

			/opt/Shared_Software/mrtrix3/bin/dwidenoise ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.nii.gz  ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise.nii.gz
			/opt/Shared_Software/mrtrix3/bin/mrdegibbs ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise.nii.gz ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs.nii.gz

			# Trova indici con bval = 0
			indices_real=($(awk '$1 == 0 {print NR}' ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvals))
			i=1
			for index in "${indices_real[@]}"
			do
	  			echo "Indice della riga con valore 0: $index"
	  			echo ${i}
	  			
	  			#estrai volume con bval = 0
				fslroi ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs.nii.gz ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs_1vol_${i}.nii.gz $((${index}-1)) 1 
				i=$((i + 1))
			done

			# Costruisci il comando fslmaths dinamicamente
			fslmaths_command="fslmaths"

			# Aggiungi ciascun file DTI al comando fslmaths
			for ((j=1; j<i; j++))
			do
	  			fslmaths_command+=" ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs_1vol_${j}.nii.gz -add"
			done

			# Rimuove l'ultima occorrenza di '-add'
			fslmaths_command="${fslmaths_command% -add}"

			# Calcola la divisione dinamica per il numero di file DTI
			fslmaths_command+=" -div $((${i}-1)) ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs_MEDIA_vols.nii.gz"

			# Esegui il comando fslmaths
			echo "Eseguendo comando: $fslmaths_command"
			eval $fslmaths_command

			#merge fra le due medie 
			fslmerge -t  ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_MEDIATED.nii.gz  ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs_MEDIA_vols.nii.gz ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_FORCORR_denoise_degibbs_MEDIA_vols.nii.gz 


			/opt/Shared_Software/mrtrix3/bin/dwifslpreproc  ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs.nii.gz ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs_PREPROC.nii.gz -rpe_pair -fslgrad  "${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvecs" "${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvals" -se_epi ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_MEDIATED.nii.gz -pe_dir AP -eddy_options " --slm=linear"

			/opt/Shared_Software/mrtrix3/bin/dwibiascorrect ants ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs_PREPROC.nii.gz ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs_PREPROC_BIASC.nii.gz -fslgrad "${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvecs" "${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvals"

		fi  

	
	fi

elif [ ${CORR} == '0' ]; then

	echo "\n\nNon eseguo TOPUP\n\n"

	#CORREZIONE EDDY PER LA DWI
	if [ ${JSON} == '1' ]; then 
		echo "\nNON Eseguo TOPUP con PE_dir presa dal .json\n"

		/opt/Shared_Software/mrtrix3/bin/dwidenoise ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.nii.gz  ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise.nii.gz
		/opt/Shared_Software/mrtrix3/bin/mrdegibbs ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise.nii.gz ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs.nii.gz

		/opt/Shared_Software/mrtrix3/bin/dwifslpreproc ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs.nii.gz ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs_PREPROC.nii.gz -rpe_none -fslgrad "${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvecs" "${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvals" -pe_dir ${PE} -eddy_options " --slm=linear"

		/opt/Shared_Software/mrtrix3/bin/dwibiascorrect ants ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs_PREPROC.nii.gz ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs_PREPROC_BIASC.nii.gz -fslgrad "${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvecs" "${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvals"

	elif [ ${JSON} == '0' ]; then
		echo "\nNON Eseguo TOPUP con PE_dir imposta\n"

		/opt/Shared_Software/mrtrix3/bin/dwidenoise ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.nii.gz  ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise.nii.gz
		/opt/Shared_Software/mrtrix3/bin/mrdegibbs ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise.nii.gz ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs.nii.gz

		/opt/Shared_Software/mrtrix3/bin/dwifslpreproc ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs.nii.gz ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs_PREPROC.nii.gz -rpe_none -fslgrad "${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvecs" "${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvals" -pe_dir AP -eddy_options " --slm=linear"

		/opt/Shared_Software/mrtrix3/bin/dwibiascorrect ants ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs_PREPROC.nii.gz ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs_PREPROC_BIASC.nii.gz -fslgrad "${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvecs" "${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}.bvals"
	fi  
fi




#REGISTRAZIONE MNI
flirt -dof 12 -in  ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs_PREPROC_BIASC.nii.gz -ref /usr/local/fsl/data/atlases/JHU/JHU-ICBM-DWI-2mm.nii.gz -omat ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_omat_${DRAMA_BPD_ID} -out ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs_PREPROC_BIASC_1vol.nii.gz

flirt -dof 12 -in ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs_PREPROC_BIASC.nii.gz -ref /usr/local/fsl/data/atlases/JHU/JHU-ICBM-DWI-2mm.nii.gz -applyxfm -init ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_omat_${DRAMA_BPD_ID} -out ${NP_PATH}/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/DTI/DTI_${DRAMA_BPD_ID}_denoise_degibbs_PREPROC_BIASC_REG.nii.gz

echo "\n\n########## END OF DTI PREPROCESSING ##########\n\n"

##################### PROCESSING - TRACULA #######################
mkdir ${WD}/PROCESSED_SUBJ/TRACULA
mkdir ${WD}/PROCESSED_SUBJ/TRACULA/${DRAMA_BPD_ID}
chmod 777 ${WD}/PROCESSED_SUBJ/TRACULA
chmod 777 ${WD}/PROCESSED_SUBJ/TRACULA/${DRAMA_BPD_ID}

cp /opt/Shared_Software/FreeSurfer-7.3.2/bin/dmrirc.example ${WD}/PROCESSED_SUBJ/TRACULA/${DRAMA_BPD_ID}/dmrirc_${DRAMA_BPD_ID}
chmod 777 ${WD}/PROCESSED_SUBJ/TRACULA/${DRAMA_BPD_ID}/dmrirc_${DRAMA_BPD_ID}


  sed -i "s=/path/to/recons/of/ducks=${WD}/PROCESSED_SUBJ/FREESURFER=g" ${WD}/PROCESSED_SUBJ/TRACULA/${DRAMA_BPD_ID}/dmrirc_${DRAMA_BPD_ID}

  sed -i "s=/path/to/tracts/of/ducks=${WD}/PROCESSED_SUBJ/TRACULA=g" ${WD}/PROCESSED_SUBJ/TRACULA/${DRAMA_BPD_ID}/dmrirc_${DRAMA_BPD_ID}
  
 # sed -i "s/set subjlist = ()/set subjlist = ()/g" ${WD}/PROCESSED_SUBJ/TRACULA/${DRAMA_BPD_ID}/dmrirc_${DRAMA_BPD_ID}
  
  sed -i "s=(huey dewey louie)=(${DRAMA_BPD_ID})=g" ${WD}/PROCESSED_SUBJ/TRACULA/${DRAMA_BPD_ID}/dmrirc_${DRAMA_BPD_ID}
  sed -i "s=/path/to/dicoms/of/ducks=${WD}/DTI=g" ${WD}/PROCESSED_SUBJ/TRACULA/${DRAMA_BPD_ID}/dmrirc_${DRAMA_BPD_ID}
 sed -i "s/dob0 = 1/dob0 = 0/g" ${WD}/PROCESSED_SUBJ/TRACULA/${DRAMA_BPD_ID}/dmrirc_${DRAMA_BPD_ID}
  
 sed -i "s/set doeddy = 2/set doeddy = 0/g" ${WD}/PROCESSED_SUBJ/TRACULA/${DRAMA_BPD_ID}/dmrirc_${DRAMA_BPD_ID}
 sed -i "s/set dTE = 0.0025/#set dTE = 0.0025/g" ${WD}/PROCESSED_SUBJ/TRACULA/${DRAMA_BPD_ID}/dmrirc_${DRAMA_BPD_ID}
  
  sed -i "s=(huey/day1/XXX-1.dcm dewey/day1/XXX-1.dcm louie/day2/XXX-1.dcm)=(DTI_${DRAMA_BPD_ID}_denoise_degibbs_PREPROC_BIASC_REG.nii.gz)=g" ${WD}/PROCESSED_SUBJ/TRACULA/${DRAMA_BPD_ID}/dmrirc_${DRAMA_BPD_ID}
 sed -i 's/set runlist = (1 3)/set runlist = (1)/g' ${WD}/PROCESSED_SUBJ/TRACULA/${DRAMA_BPD_ID}/dmrirc_${DRAMA_BPD_ID}
  
  sed -i "65s#.*#set bveclist = (DTI_${DRAMA_BPD_ID}.bvecs)#g" ${WD}/PROCESSED_SUBJ/TRACULA/${DRAMA_BPD_ID}/dmrirc_${DRAMA_BPD_ID}
 sed -i '66s#.*##g' ${WD}/PROCESSED_SUBJ/TRACULA/${DRAMA_BPD_ID}/dmrirc_${DRAMA_BPD_ID}
  
 sed -i '67s#.*##g' ${WD}/PROCESSED_SUBJ/TRACULA/${DRAMA_BPD_ID}/dmrirc_${DRAMA_BPD_ID}
  
  sed -i "74s#.*#set bvallist = (DTI_${DRAMA_BPD_ID}.bvals)#g" ${WD}/PROCESSED_SUBJ/TRACULA/${DRAMA_BPD_ID}/dmrirc_${DRAMA_BPD_ID}
 sed -i '75s#.*##g' ${WD}/PROCESSED_SUBJ/TRACULA/${DRAMA_BPD_ID}/dmrirc_${DRAMA_BPD_ID}
  
 sed -i '76s#.*##g' ${WD}/PROCESSED_SUBJ/TRACULA/${DRAMA_BPD_ID}/dmrirc_${DRAMA_BPD_ID}
  
 sed -i 's!set b0mlist = (huey/fmag/XXX-1.dcm dewey/fmag/XXX-1.dcm louie/fmag/XXX-1.dcm)!#set b0mlist = (huey/fmag/XXX-1.dcm dewey/fmag/XXX-1.dcm louie/fmag/XXX-1.dcm)!g' ${WD}/PROCESSED_SUBJ/TRACULA/${DRAMA_BPD_ID}/dmrirc_${DRAMA_BPD_ID}
  
 sed -i 's!set b0plist = (huey/fphas/XXX-1.dcm dewey/fphas/XXX-1.dcm louie/fphas/XXX-1.dcm)!#set b0plist = (huey/fphas/XXX-1.dcm dewey/fphas/XXX-1.dcm louie/fphas/XXX-1.dcm)!g' ${WD}/PROCESSED_SUBJ/TRACULA/${DRAMA_BPD_ID}/dmrirc_${DRAMA_BPD_ID}
  
 sed -i 's!set pathlist = ( lh.uf rh.uf cc.rostrum )!set pathlist = ( acomm cc.bodyc cc.bodyp cc.bodypf cc.bodypm cc.bodyt cc.genu cc.rostrum cc.splenium mcp lh.af rh.af lh.ar rh.ar lh.atr rh.atr lh.cbd rh.cbd lh.cbv rh.cbv rh.cst lh.emc rh.emc lh.cst lh.fat rh.fat lh.fx rh.fx lh.ilf rh.ilf lh.mlf rh.mlf lh.or rh.or lh.slf1 rh.slf1 lh.slf2 rh.slf2 lh.slf3 rh.slf3 lh.uf rh.uf )!g' ${WD}/PROCESSED_SUBJ/TRACULA/${DRAMA_BPD_ID}/dmrirc_${DRAMA_BPD_ID}
  
 sed -i 's!set ncpts = ( 7 7 5 )!set ncpts = ( 7 8 9 7 7 12 5 5 7 7 9 9 4 4 4 4 7 7 5 5 7 7 7 7 5 5 9 9 8 8 7 7 5 5 7 7 8 8 6 6 7 7 )!g' ${WD}/PROCESSED_SUBJ/TRACULA/${DRAMA_BPD_ID}/dmrirc_${DRAMA_BPD_ID}
  
 sed -i 's!set intertrg = /path/to/a/duck/template.nii.gz!#set intertrg = /path/to/a/duck/template.nii.gz!g' ${WD}/PROCESSED_SUBJ/TRACULA/${DRAMA_BPD_ID}/dmrirc_${DRAMA_BPD_ID}
  
 sed -i 's!set intradof = 6!set intradof = 12!g' ${WD}/PROCESSED_SUBJ/TRACULA/${DRAMA_BPD_ID}/dmrirc_${DRAMA_BPD_ID}
  
 sed -i 's!set thrbet = 0.5!set thrbet = 0.3!g' ${WD}/PROCESSED_SUBJ/TRACULA/${DRAMA_BPD_ID}/dmrirc_${DRAMA_BPD_ID}
  
 sed -i 's!set echospacing = 0.7!#set echospacing = 0.7!g' ${WD}/PROCESSED_SUBJ/TRACULA/${DRAMA_BPD_ID}/dmrirc_${DRAMA_BPD_ID}


trac-all -prep -c ${WD}/PROCESSED_SUBJ/TRACULA/${DRAMA_BPD_ID}/dmrirc_${DRAMA_BPD_ID}
trac-all -bedp -c ${WD}/PROCESSED_SUBJ/TRACULA/${DRAMA_BPD_ID}/dmrirc_${DRAMA_BPD_ID}
trac-all -path -c ${WD}/PROCESSED_SUBJ/TRACULA/${DRAMA_BPD_ID}/dmrirc_${DRAMA_BPD_ID}

# estrazione risultati 
mkdir ${WD}/PROCESSED_SUBJ/TRACULA/${DRAMA_BPD_ID}/Output_txt

for aTract in acomm_avg16_syn_bbr  cc.bodyc_avg16_syn_bbr cc.bodyp_avg16_syn_bbr cc.bodypf_avg16_syn_bbr cc.bodypm_avg16_syn_bbr cc.bodyt_avg16_syn_bbr cc.genu_avg16_syn_bbr cc.rostrum_avg16_syn_bbr cc.splenium_avg16_syn_bbr   mcp_avg16_syn_bbr lh.af_avg16_syn_bbr rh.af_avg16_syn_bbr lh.ar_avg16_syn_bbr rh.ar_avg16_syn_bbr lh.atr_avg16_syn_bbr rh.atr_avg16_syn_bbr     lh.cbd_avg16_syn_bbr rh.cbd_avg16_syn_bbr lh.cbv_avg16_syn_bbr rh.cbv_avg16_syn_bbr rh.cst_avg16_syn_bbr lh.emc_avg16_syn_bbr rh.emc_avg16_syn_bbr lh.cst_avg16_syn_bbr     lh.fat_avg16_syn_bbr rh.fat_avg16_syn_bbr lh.fx_avg16_syn_bbr rh.fx_avg16_syn_bbr lh.ilf_avg16_syn_bbr rh.ilf_avg16_syn_bbr lh.mlf_avg16_syn_bbr rh.mlf_avg16_syn_bbr     lh.or_avg16_syn_bbr rh.or_avg16_syn_bbr lh.slf1_avg16_syn_bbr rh.slf1_avg16_syn_bbr lh.slf2_avg16_syn_bbr rh.slf2_avg16_syn_bbr lh.slf3_avg16_syn_bbr rh.slf3_avg16_syn_bbr lh.uf_avg16_syn_bbr rh.uf_avg16_syn_bbr; do echo $aTract; tractstats2table --inputs ${WD}/PROCESSED_SUBJ/TRACULA/${DRAMA_BPD_ID}/dpath/${aTract}/pathstats.overall.txt  --overall  --only-measures FA_Avg_Weight MD_Avg_Weight -t ${WD}/PROCESSED_SUBJ/TRACULA/${DRAMA_BPD_ID}/Output_txt/${aTract}_${DRAMA_BPD_ID}.txt -v; sed 's/ \+/,/g' ${WD}/PROCESSED_SUBJ/TRACULA/${DRAMA_BPD_ID}/Output_txt/${aTract}_${DRAMA_BPD_ID}.txt > ${WD}/PROCESSED_SUBJ/TRACULA/${DRAMA_BPD_ID}/Output_txt/${aTract}_${DRAMA_BPD_ID}.csv; sed -e 's#MD_Avg_Weight#'MD_Avg_Weight_${aTract}'#g' -e 's#FA_Avg_Weight#'FA_Avg_Weight_${aTract}'#g' ${WD}/PROCESSED_SUBJ/TRACULA/${DRAMA_BPD_ID}/Output_txt/${aTract}_${DRAMA_BPD_ID}.csv > ${WD}/PROCESSED_SUBJ/TRACULA/${DRAMA_BPD_ID}/Output_txt/${aTract}_${DRAMA_BPD_ID}_name.csv; cut -f2- ${WD}/PROCESSED_SUBJ/TRACULA/${DRAMA_BPD_ID}/Output_txt/${aTract}_${DRAMA_BPD_ID}_name.csv >  ${WD}/PROCESSED_SUBJ/TRACULA/${DRAMA_BPD_ID}/Output_txt/${aTract}_${DRAMA_BPD_ID}_name_finale.csv; rm ${WD}/PROCESSED_SUBJ/TRACULA/${DRAMA_BPD_ID}/Output_txt/${aTract}_${DRAMA_BPD_ID}_name.csv; done

paste -d, ${WD}/PROCESSED_SUBJ/TRACULA/${DRAMA_BPD_ID}/Output_txt/*name_finale.csv > ${RESULTS_DIR}/${DRAMA_BPD_ID}_Output_TRACULA.csv

#awk -F '\t' -v OFS=',' '{print $1}' "${RESULTS_DIR}/${DRAMA_BPD_ID}_Output_TRACULA.csv" > ${WD}/RESULTS_DIR/${DRAMA_BPD_ID}_Output_TRACULA_features_selected.csv
elif [ "${FileName_dti}" == "NULL" ]; then
echo "DTI non presente"
fi

paste -d, ${RESULTS_DIR}/${DRAMA_BPD_ID}_Output_FS_TOTALE.csv  ${RESULTS_DIR}/${DRAMA_BPD_ID}_Output_TRACULA.csv > ${RESULTS_DIR}/${DRAMA_BPD_ID}_Output_TOTALE_IMAGING.csv

rm ${WD}/pipelinesrunning
touch ${WD}/pipelinesfinished


# Estrazione Features significative

sed 's/\t/,/g' ${RESULTS_DIR}/${DRAMA_BPD_ID}_Output_TOTALE_IMAGING.csv > ${RESULTS_DIR}/${DRAMA_BPD_ID}_Output_TOTALE_IMAGING_modified.csv
csvcut -c MD_Avg_Weight_mcp_avg16_syn_bbr,Left-Inf-Lat-Vent,lh_caudalmiddlefrontal_thickness,rh_caudalmiddlefrontal_thickness,lh_hippocampal-fissure,CC_Mid_Anterior ${RESULTS_DIR}/${DRAMA_BPD_ID}_Output_TOTALE_IMAGING_modified.csv > ${RESULTS_DIR}/${DRAMA_BPD_ID}_Output_IMAGING_SELECTED_FEATURES.csv

LeftInfLatVent=$(csvcut -c Left-Inf-Lat-Vent "${RESULTS_DIR}/${DRAMA_BPD_ID}_Output_IMAGING_SELECTED_FEATURES.csv" | tail -n +2)
CCMidAnterior=$(csvcut -c CC_Mid_Anterior "${RESULTS_DIR}/${DRAMA_BPD_ID}_Output_IMAGING_SELECTED_FEATURES.csv" | tail -n +2)
lhcaudalmiddlefrontalthickness=$(csvcut -c lh_caudalmiddlefrontal_thickness "${RESULTS_DIR}/${DRAMA_BPD_ID}_Output_IMAGING_SELECTED_FEATURES.csv" | tail -n +2)
rhcaudalmiddlefrontalthickness=$(csvcut -c rh_caudalmiddlefrontal_thickness "${RESULTS_DIR}/${DRAMA_BPD_ID}_Output_IMAGING_SELECTED_FEATURES.csv" | tail -n +2)
lhhippocampalfissure=$(csvcut -c lh_hippocampal-fissure "${RESULTS_DIR}/${DRAMA_BPD_ID}_Output_IMAGING_SELECTED_FEATURES.csv" | tail -n +2)
MDAvgWeightmcpavg16synbbr=$(csvcut -c MD_Avg_Weight_mcp_avg16_syn_bbr "${RESULTS_DIR}/${DRAMA_BPD_ID}_Output_IMAGING_SELECTED_FEATURES.csv" | tail -n +2)


rm ${RESULTS_DIR}/${DRAMA_BPD_ID}_Output_TOTALE_IMAGING_modified.csv

#################### Per Creazione Report ####################

# masks for report
cd ${WD}/PROCESSED_SUBJ/FREESURFER/${DRAMA_BPD_ID}
cp ${NP_PATH}/commons/DRAMA_BPD/IMG_ASEG.r ${WD}/PROCESSED_SUBJ/FREESURFER/${DRAMA_BPD_ID}
mkdir masks

# ASEG #
#Left-Inf-Lat_Vent
mri_binarize --i mri/aseg.mgz --o masks/Left-Inf-Lat-Vent.nii.gz --match 5
cog_LILV=$(fslstats ${WD}/PROCESSED_SUBJ/FREESURFER/${DRAMA_BPD_ID}/masks/Left-Inf-Lat-Vent.nii.gz -C)
cog_LILV=$(echo ${cog_LILV// /,})
cog_LILV=$(echo $cog_LILV | sed 's/\(.*\),/\1 /')
sed -i "s/@@cogLILV@@/${cog_LILV}/g" IMG_ASEG.r

#CC_Mid_Anterior
mri_binarize --i mri/aseg.mgz --o masks/CC_Mid_Anterior.nii.gz --match 254
cog_CCMA=$(fslstats ${WD}/PROCESSED_SUBJ/FREESURFER/${DRAMA_BPD_ID}/masks/CC_Mid_Anterior.nii.gz -C)
cog_CCMA=$(echo ${cog_CCMA// /,})
cog_CCMA=$(echo $cog_CCMA | sed 's/\(.*\),/\1 /')
sed -i "s/@@cogCCMA@@/${cog_CCMA}/g" IMG_ASEG.r

#lh_hippocampal_fissure
mri_binarize --i mri/lh.hippoAmygLabels-T1.v22.FSvoxelSpace.mgz --o masks/lh_hippocampal_fissure.nii.gz --match 215
cog_LHF=$(fslstats ${WD}/PROCESSED_SUBJ/FREESURFER/${DRAMA_BPD_ID}/masks/lh_hippocampal_fissure.nii.gz -C)
cog_LHF=$(echo ${cog_LHF// /,})
cog_LHF=$(echo $cog_LHF | sed 's/\(.*\),/\1 /')
sed -i "s/@@cogLHF@@/${cog_LHF}/g" IMG_ASEG.r

chmod +x IMG_ASEG.r
./IMG_ASEG.r

convert Left-Inf-Lat-Vent.png -crop +35+40 -crop -265-260 Left-Inf-Lat-Vent_1.png
convert CC_Mid_Anterior.png -crop +35+40 -crop -265-260 CC_Mid_Anterior_1.png
convert lh_hippocampal_fissure.png -crop +35+40 -crop -265-260 lh_hippocampal_fissure_1.png


convert Left-Inf-Lat-Vent.png -crop +260+40 -crop -40-260 -rotate 90 Left-Inf-Lat-Vent_2.png
convert CC_Mid_Anterior.png -crop +260+40 -crop -40-260 -rotate 90 CC_Mid_Anterior_2.png
convert lh_hippocampal_fissure.png -crop +260+40 -crop -40-260 -rotate 90 lh_hippocampal_fissure_2.png


convert Left-Inf-Lat-Vent.png -crop +30+270 -crop -270-30 -rotate 180 Left-Inf-Lat-Vent_3.png
convert CC_Mid_Anterior.png -crop +30+270 -crop -270-30 -rotate 180 CC_Mid_Anterior_3.png
convert lh_hippocampal_fissure.png -crop +30+270 -crop -270-30 -rotate 180 lh_hippocampal_fissure_3.png

# APARC #
cp ${NP_PATH}/commons/DRAMA_BPD/IMG_APARC.r ${WD}/PROCESSED_SUBJ/FREESURFER/${DRAMA_BPD_ID}

#lh_caudalmiddlefrontal_thickness
mri_binarize --i mri/aparc+aseg.mgz --o masks/lh_caudalmiddlefrontal_thickness.nii.gz --match 1003
cog_LCMFT=$(fslstats ${WD}/PROCESSED_SUBJ/FREESURFER/${DRAMA_BPD_ID}/masks/lh_caudalmiddlefrontal_thickness.nii.gz -C)
cog_LCMFT=$(echo ${cog_LCMFT// /,})
cog_LCMFT=$(echo $cog_LCMFT | sed 's/\(.*\),/\1 /')
sed -i "s/@@cogLCMFT@@/${cog_LCMFT}/g" IMG_APARC.r

#rh_caudalmiddlefrontal_thickness
mri_binarize --i mri/aparc+aseg.mgz --o masks/rh_caudalmiddlefrontal_thickness.nii.gz --match 2003
cog_RCMFT=$(fslstats ${WD}/PROCESSED_SUBJ/FREESURFER/${DRAMA_BPD_ID}/masks/rh_caudalmiddlefrontal_thickness.nii.gz -C)
cog_RCMFT=$(echo ${cog_RCMFT// /,})
cog_RCMFT=$(echo $cog_RCMFT | sed 's/\(.*\),/\1 /')
sed -i "s/@@cogRCMFT@@/${cog_RCMFT}/g" IMG_APARC.r

chmod +x IMG_APARC.r
./IMG_APARC.r

convert lh_caudalmiddlefrontal_thickness.png -crop +35+40 -crop -265-260 lh_caudalmiddlefrontal_thickness_1.png
convert lh_caudalmiddlefrontal_thickness.png -crop +260+40 -crop -40-260 -rotate 90 lh_caudalmiddlefrontal_thickness_2.png
convert lh_caudalmiddlefrontal_thickness.png -crop +30+270 -crop -270-30 -rotate 180 lh_caudalmiddlefrontal_thickness_3.png

convert rh_caudalmiddlefrontal_thickness.png -crop +35+40 -crop -265-260 rh_caudalmiddlefrontal_thickness_1.png
convert rh_caudalmiddlefrontal_thickness.png -crop +260+40 -crop -40-260 -rotate 90 rh_caudalmiddlefrontal_thickness_2.png
convert rh_caudalmiddlefrontal_thickness.png -crop +30+270 -crop -270-30 -rotate 180 rh_caudalmiddlefrontal_thickness_3.png

# TRACULA #
cd ${WD}/PROCESSED_SUBJ/TRACULA/${DRAMA_BPD_ID}/dpath
cp ${NP_PATH}/commons/DRAMA_BPD/IMG_TRACULA.r ${WD}/PROCESSED_SUBJ/TRACULA/${DRAMA_BPD_ID}/dpath

#MD_Avg_Weight_mcp_avg16_syn_bbr
cogcMDMCP=$(fslstats mcp_avg16_syn_bbr/path.pd.nii.gz -C)
cogcMDMCP=$(echo ${cogcMDMCP// /,})
cogcMDMCP=$(echo $cogcMDMCP | sed 's/\(.*\),/\1 /')
sed -i "s/@@cogcMDMCP@@/${cogcMDMCP}/g" IMG_TRACULA.r

chmod +x IMG_TRACULA.r
./IMG_TRACULA.r

convert MDMCP.png -crop +35+40 -crop -265-260 MDMCP_1.png
convert MDMCP.png -crop +260+40 -crop -40-260 MDMCP_2.png
convert MDMCP.png -crop +30+270 -crop -270-30 -rotate 180 MDMCP_3.png


#################### NEUROHARMONIZE ###################

mkdir $WD/PROCESSED_SUBJ/NEUROHARMONIZE
chmod 777 -R $WD/PROCESSED_SUBJ/NEUROHARMONIZE
cp ${RESULTS_DIR}/${DRAMA_BPD_ID}_Output_IMAGING_SELECTED_FEATURES.csv $WD/PROCESSED_SUBJ/NEUROHARMONIZE/${DRAMA_BPD_ID}_Output_IMAGING_SELECTED_FEATURES.csv
cp ${NP_PATH}/commons/DRAMA_BPD/IMAGING_RR_QC_LIGHT_VERSION_for_Suicide_Squad_per_armonizzazione_TOT_only_CLM.csv $WD/PROCESSED_SUBJ/NEUROHARMONIZE/IMAGING_RR_QC_LIGHT_VERSION_for_Suicide_Squad_per_armonizzazione_TOT_only_CLM.csv
cp ${NP_PATH}/commons/DRAMA_BPD/covariates_TOT_only_CLM.csv $WD/PROCESSED_SUBJ/NEUROHARMONIZE/covariates_TOT_only_CLM.csv
cp ${NP_PATH}/commons/DRAMA_BPD/MY_MODEL_* $WD/PROCESSED_SUBJ/NEUROHARMONIZE/
cp ${NP_PATH}/commons/DRAMA_BPD/Harmonization_DRAMA_BPD.py $WD/PROCESSED_SUBJ/NEUROHARMONIZE/Harmonization_DRAMA_BPD.py

cd $WD/PROCESSED_SUBJ/NEUROHARMONIZE

virtualenv HARM_VE
source HARM_VE/bin/activate
pip install -r requirements_HARM.txt

python3 Harmonization_DRAMA_BPD.py --ID ${DRAMA_BPD_ID}

deactivate

rm $WD/PROCESSED_SUBJ/NEUROHARMONIZE/IMAGING_RR_QC_LIGHT_VERSION_for_Suicide_Squad_per_armonizzazione_TOT_only_CLM.csv
rm $WD/PROCESSED_SUBJ/NEUROHARMONIZE/covariates_TOT_only_CLM.csv
rm $WD/PROCESSED_SUBJ/NEUROHARMONIZE/MY_MODEL_*


################### CLASSIFIER ########################

mkdir $WD/PROCESSED_SUBJ/CLASSIFIER
chmod 777 -R $WD/PROCESSED_SUBJ/CLASSIFIER
cd $WD/PROCESSED_SUBJ/CLASSIFIER

mv $WD/PROCESSED_SUBJ/NEUROHARMONIZE/New_subject_harmonized.csv $WD/PROCESSED_SUBJ/NEUROHARMONIZE/${DRAMA_BPD_ID}_harmonized.csv
cp $WD/PROCESSED_SUBJ/NEUROHARMONIZE/${DRAMA_BPD_ID}_harmonized.csv $WD/PROCESSED_SUBJ/CLASSIFIER/Imaging_UC.csv

#Trasforma maiuscole in minuscole per classificatore

file_csv="Imaging_UC.csv"
output_csv="Imaging.csv"
cat "$file_csv" | tr '[:upper:]' '[:lower:]' > "$output_csv"

cp ${NP_PATH}/commons/DRAMA_BPD/DRAMA-BPD.py $WD/PROCESSED_SUBJ/CLASSIFIER
cp ${NP_PATH}/commons/DRAMA_BPD/DRAMA-BPD.pkl $WD/PROCESSED_SUBJ/CLASSIFIER
cp ${NP_PATH}/commons/DRAMA_BPD/DRAMA-BPD_SCALER.pkl $WD/PROCESSED_SUBJ/CLASSIFIER
cp ${NP_PATH}/commons/DRAMA_BPD/requirements.txt $WD/PROCESSED_SUBJ/CLASSIFIER


# Aggiungi variabili cliniche
var1="farmnp_t0"
var2="scl90_tot_t0"
clin1=${farmnp_t0}
clin2=${scl90_tot_t0}

clin_csv="clin.csv"

# Crea il file CSV e aggiungi le intestazioni (i nomi delle variabili)
echo "$var1,$var2" > "$clin_csv"

# Aggiungi i valori delle variabili
echo "$clin1,$clin2" >> "$clin_csv"

# Merge csv imaging + clinico
file1="Imaging.csv"
file2="clin.csv"
Input_csv="Input.csv"

paste -d',' "$file1" "$file2" > "$Input_csv"

# Crea e lavora nell'ambiente virtuale

virtualenv DRAMA_BPD_VE
source DRAMA_BPD_VE/bin/activate
pip install -r requirements.txt

python3 DRAMA-BPD.py

deactivate


# Nome del file CSV
output_csv="Output.csv"

# Leggi il primo valore dalla prima riga e prima colonna
ts=$(awk -F',' 'NR==1 {print $1}' "$output_csv")

AUC=75

if [ ${ts} == '1.0' ]; then 
	
	classe=SA
	echo "Soggetto ${DRAMA_BPD_ID} classificato come: ${classe} con AUC = ${AUC}%"
elif [ ${ts} == '0.0' ]; then
	
	classe=NA
	echo "Soggetto ${DRAMA_BPD_ID} classificato come: ${classe} con AUC = ${AUC}%"
	AUC=$((-${AUC}))
	
fi


##################### LATEX REPORT ####################

cd $WD/PROCESSED_SUBJ/CLASSIFIER
cp ${NP_PATH}/commons/DRAMA_BPD/CREATE_PLOT.r $WD/PROCESSED_SUBJ/CLASSIFIER

# Plotting results with R
echo "Plotting results with R"
sed -i "s/@@AUC@@/${AUC}/g" CREATE_PLOT.r

chmod +x CREATE_PLOT.r
./CREATE_PLOT.r

cp ${NP_PATH}/commons/DRAMA_BPD/template.tex ${WD}/REPORT_DRAMA_DBP.tex

cd ${WD}

sed -i "s#@@ID@@#$DRAMA_BPD_ID#" REPORT_DRAMA_DBP.tex
sed -i "s#@@ID@@#$DRAMA_BPD_ID#" REPORT_DRAMA_DBP.tex
sed -i "s#@@ID@@#$DRAMA_BPD_ID#" REPORT_DRAMA_DBP.tex
sed -i "s#@@LeftInfLatVent@@#$LeftInfLatVent#" REPORT_DRAMA_DBP.tex
sed -i "s#@@CCMidAnterior@@#$CCMidAnterior#" REPORT_DRAMA_DBP.tex
sed -i "s#@@lhhippocampalfissure@@#$lhhippocampalfissure#" REPORT_DRAMA_DBP.tex
sed -i "s#@@lhcaudalmiddlefrontalthickness@@#$lhcaudalmiddlefrontalthickness#" REPORT_DRAMA_DBP.tex
sed -i "s#@@rhcaudalmiddlefrontalthickness@@#$rhcaudalmiddlefrontalthickness#" REPORT_DRAMA_DBP.tex
sed -i "s#@@MDAvgWeightmcpavg16synbbr@@#$MDAvgWeightmcpavg16synbbr#" REPORT_DRAMA_DBP.tex
sed -i "s#@@classe@@#${classe}#" REPORT_DRAMA_DBP.tex

pdflatex REPORT_DRAMA_DBP.tex
mv REPORT_DRAMA_DBP.pdf ${DRAMA_BPD_ID}.pdf

################### SENDING REPORT ########################

cd $WD
mkdir Output


# Create file for download output
mv Output /media/NewPsy/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/Output_ng2
#FS
cp $WD/PROCESSED_SUBJ/FREESURFER/${DRAMA_BPD_ID}/mri/aparc+aseg.mgz /media/NewPsy/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/Output_ng2
cp $WD/PROCESSED_SUBJ/FREESURFER/${DRAMA_BPD_ID}/mri/orig.mgz /media/NewPsy/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/Output_ng2
#TRACULA
cp $WD/PROCESSED_SUBJ/TRACULA/${DRAMA_BPD_ID}/dmri/dtifit_FA.nii.gz /media/NewPsy/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/Output_ng2
cp $WD/PROCESSED_SUBJ/TRACULA/${DRAMA_BPD_ID}/dpath/merged_avg16_syn_bbr.mgz /media/NewPsy/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/Output_ng2
#RESULTS
cp ${RESULTS_DIR}/${DRAMA_BPD_ID}_Output_TOTALE_IMAGING.csv /media/NewPsy/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/Output_ng2
cp ${RESULTS_DIR}/${DRAMA_BPD_ID}_Output_IMAGING_SELECTED_FEATURES.csv /media/NewPsy/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/Output_ng2
cp $WD/PROCESSED_SUBJ/NEUROHARMONIZE/${DRAMA_BPD_ID}_harmonized.csv /media/NewPsy/sandbox/${USER_NAME_ID}/${job_folder}/${DRAMA_BPD_ID}/Output_ng2


zip -r Output_ng2.zip Output_ng2
rm -rf Output_ng2


# Sending report through Mail
cp /media/NewPsy/script_newpsy/mail.sh ./mail.sh
chmod +x ./mail.sh

./mail.sh "DRAMA_BPD" "${WD}/${DRAMA_BPD_ID}.pdf" "${user_email}" "${WD}"





















