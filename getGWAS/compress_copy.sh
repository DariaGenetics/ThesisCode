#!/usr/bin/env bash
#
#SBATCH -J bam_qc # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 5:00:00 ###
#SBATCH --mem 2
#SBATCH -o /scratch/aob2x/dest/slurmOutput/bam_qc.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/bam_qc.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard

# sbatch ~/ThesisCode/getGWAS/compress_copy.sh
# sacct -j 46191468

tar -czvf /scratch/aob2x/gwas.test.tar.gz \
/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/*/ChillComaRecoveryTime_standard_female*/*txt \
/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/*/ChillComaRecoveryTime_standard_male*/*txt \
/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/*/ClimbingHeight_ParaquatInsecticideExposure_M*/*txt \
/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/*/ClimbingHeight_standard_M*/*txt \
/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/*/DayBoutNumber_standard_Female*/*txt \
/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/*/DaySleep_standard_Female*/*txt \
/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/*/FoodIntake_standard_female*/*txt \
/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/*/FreeGlycerolLevels_HighGlucoseDiet_male*/*txt \
/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/*/FreeGlycerolLevels_LowGlucoseDiet_male*/*txt \
/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/*/FreeGlycerolLevels_PooledDiets_male*/*txt \
/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/*/LarvaeSurvival_10Î¼g-mLChlorantraniliproleOn10mLFoodMedia_LarvaeMixed*/*txt \
/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/*/LarvaeSurvival_20Î¼g-mLChlorantraniliproleOn10mLFoodMedia_LarvaeMixed*/*txt \
/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/*/LarvaeSurvival_40Î¼g-mLChlorantraniliproleOn10mLFoodMedia_LarvaeMixed*/*txt \
/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/*/MeanElutionTime_RepeatedEthanolExposures_female*/*txt \
/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/*/MeanElutionTime_RepeatedEthanolExposures_male*/*txt \
/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/*/MinimumSurvivableTemperature_TwentyNineDegreesC_M*/*txt \
/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/*/MinimumSurvivableTemperature_TwentySixDegreesC_M*/*txt \
/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/*/MinimumSurvivableTemperature_TwentyThreeDegreesC_M*/*txt \
/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/*/SolubleProteinLevels_LowGlucoseDiet_male*/*txt \
/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/*/SolubleProteinLevels_PooledDiets_male*/*txt \
/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/*/Weight_HighGlucoseDiet_male*/*txt \
/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/*/Weight_LowGlucoseDiet_male*/*txt \
/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/*/Weight_PooledDiets_male*/*txt \
/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/*/1-HexanolAttraction_standard_female*/*txt \
/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/*/2-HeptanoneAttraction_standard_female*/*txt \
/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/*/2-HeptanoneAttraction_standard_male*/*txt \
/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/*/CentroidSize_standard_F*/*txt \
/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/*/CentroidSize_standard_M*/*txt \
/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/*/HexanalAttraction_standard_female*/*txt \
/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/*/HexanalAttraction_standard_male*/*txt \
/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/*/IntensityOfResponse_standard_F*/*txt \
/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/*/InternocularDistance_standard_F*/*txt \
/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/*/InternocularDistance_standard_M*/*txt \
/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/*/InterocularDistance_standard_female*/*txt \
/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/*/InterocularDistance_standard_male*/*txt \
/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/*/LifeTimeFecundity_LowYeast_female*/*txt \
/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/*/NegativeGeotaxis_MSB_Treatment_female*/*txt \
/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/*/PupalCaseLength_Standard_Mix*/*txt \
/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/*/ThoraxLength_standard_female*/*txt \
/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/*/TriglycerideLevels_LowGlucoseDiet_male*/*txt \
/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/*/TriglycerideLevels_PooledDiets_male*/*txt \
/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/*/Weight_HighGlucoseDiet_male*/*txt \
/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/*/Weight_LowGlucoseDiet_male*/*txt \
/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/*/Weight_PooledDiets_male*/*txt \
/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/*/WingCentroidSize_standard_female*/*txt \
/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/*/WingCentroidSize_standard_male*/*txt
