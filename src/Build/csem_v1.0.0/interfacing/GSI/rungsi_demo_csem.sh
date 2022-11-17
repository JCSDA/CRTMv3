#!/bin/sh --login

#PBS -A da-cpu 
#PBS -N global_gsi
#PBS -j oe
#PBS -q batch 
#PBS -l nodes=50:ppn=6
#PBS -l walltime=00:30:00

set -x

ulimit -s unlimited
ulimit -a

. $MODULESHOME/init/sh
##module load PrgEnv-intel
##module load cray-mpich
##module load prod_util
##module load prod_envir
##module load crtm-intel/2.2.4
##module load cfp-intel-sandybridge
module list

export OMP_STACKSIZE=2G
export OMP_NUM_THREADS=4
export FORT_BUFFERED=true

# Set experiment name and analysis date
adate=2018020100

exp=ProdGSI-csem.$adate


# Set YES to use ensemble, NO=standard 3dvar
DOHYBVAR=NO

# Set YES to run 4D-EnsVar.  NO=3D-EnsVar or 3DVAR
DO4DENSVAR=NO

# Set new radiance bias correction flag
export UNCOMPRESS=gunzip

# Generate diagnostic files
GENDIAG=YES
DIAG_SUFFIX=""
CDATE=$adate
DIAG_COMPRESS=YES
COMPRESS=gzip
DIAG_TARBALL=YES
USE_CFP=NO


# Select data dump (gdas creates diagnostic files, gfs has no diagnostic files)
dumpobs=gdas


# Size of ensemble
ENS_NUM_ANAL=80


# Set path/file for gsi executable
#gsiexec=/scratch4/NCEPDEV/da/save/Emily.Liu/GSI_Demo/ProdGSI/src/global_gsi
gsiexec=/scratch4/NCEPDEV/jcsda/save/Ming.Chen/ProdGSI/bin/global_gsi

# Set the JCAP resolution which you want.
# All resolutions use LEVS=64
export JCAP=574
export JCAP_B=1534
export LEVS=64

# Set runtime and save directories
PTMP=/scratch4/NCEPDEV/stmp3
DATA=$PTMP/$LOGNAME/tmp${JCAP}_sigmap/${exp}
SAVDIR=$PTMP/$LOGNAME/out${JCAP}_sigmap/${exp}



# Specify GSI fixed field
#fixgsi=/scratch4/NCEPDEV/da/save/Emily.Liu/GSI_Demo/ProdGSI/fix
#fixcrtm=/scratch4/NCEPDEV/da/save/Michael.Lueken/nwprod/lib/crtm/2.2.3/fix_update
fixgsi=/scratch4/NCEPDEV/da/save/Emily.Liu/GSI_Demo/ProdGSI/fix
fixcrtm=/scratch4/NCEPDEV/jcsda/save/Ming.Chen/CRTM_CSEM/crtm_v2.3.0/fix 


# Set variables used in script
#   CLEAN up $DATA when finished (YES=remove, NO=leave alone)
#   ndate is a date manipulation utility
#   ncp is cp replacement, currently keep as /bin/cp

CLEAN=NO
ndate=${NDATE:-/scratch4/NCEPDEV/da/save/Michael.Lueken/nwprod/util/exec/ndate}
export wc=${wc:-/usr/bin/wc}
ncp1=/bin/cp
ncp="ln -fs"

# Given the requested resolution, set dependent resolution parameters
if [[ "$JCAP" = "1534" ]]; then
   export LONA=3072
   export LATA=1536
   export DELTIM=120
   export resol=1
elif [[ "$JCAP" = "1148" ]]; then
   export LONA=2304
   export LATA=1152
   export DELTIM=120
   export resol=1
elif [[ "$JCAP" = "878" ]]; then
   export LONA=1760
   export LATA=880
   export DELTIM=120
   export resol=1
elif [[ "$JCAP" = "574" ]]; then
   export LONA=1152
   export LATA=576
   export DELTIM=450
   export resol=1
elif [[ "$JCAP" = "382" ]]; then
   export LONA=768
   export LATA=384
   export DELTIM=180
   export resol=1
elif [[ "$JCAP" = "254" ]]; then
   export LONA=512
   export LATA=256
   export DELTIM=450
   export resol=1
elif [[ "$JCAP" = "126" ]]; then
   export LONA=384
   export LATA=190
   export DELTIM=600
   export resol=1
elif [[ "$JCAP" = "62" ]]; then
   export LONA=192
   export LATA=94
   export DELTIM=1200
   export resol=2
else
   echo "INVALID JCAP = $JCAP"
   exit
fi
export NLON=$LONA
export NLAT=$((${LATA}+2))

# Given the analysis date, compute the date from which the
# first guess comes.  Extract cycle and set prefix and suffix
# for guess and observation data files
gdate=`$ndate -06 $adate`
PDYa=`echo $adate | cut -c1-8`
cyca=`echo $adate | cut -c9-10`
PDYg=`echo $gdate | cut -c1-8`
cycg=`echo $gdate | cut -c9-10`
if [ $dumpobs = gdas ]; then
   prefix_obs=gdas.t${cyca}z
else
   prefix_obs=gfs.t${cyca}z
fi
prefix_tbc=gdas.t${cycg}z
prefix_sfc=gdas.t${cycg}z
prefix_atm=gdas.t${cyca}z
prefix_ens=gdas.t${cycg}z
suffix=tm00.bufr_d

adate0=`echo $adate | cut -c1-8`
gdate0=`echo $gdate | cut -c1-8`

dumpges=gdas
dumpatm=gdas


datobs=/scratch4/NCEPDEV/stmp3/$USER/data_sigmap/globalprod.${CDATE}
datges=/scratch4/NCEPDEV/stmp3/$USER/data_sigmap/globalprod.${CDATE}
datatm=$datges
datens=$datges

# Set up $DATA
rm -rf $DATA
mkdir -p $DATA
cd $DATA
rm -rf core*

# Make gsi namelist

# CO2 namelist and file decisions
ICO2=${ICO2:-0}
if [ $ICO2 -gt 0 ] ; then
#  Copy co2 files to $DATA
   co2dir=${CO2DIR:-$fixgsi}
   yyyy=$(echo ${CDATE:-$adate}|cut -c1-4)
   rm ./global_co2_data.txt
   co2=$co2dir/global_co2.gcmscl_$yyyy.txt
   if [ -s $co2 ] ; then
      $ncp1 $co2 ./global_co2_data.txt
   fi
   if [ ! -s ./global_co2_data.txt ] ; then
      echo "\./global_co2_data.txt" not created
##    exit 1
   fi
fi

SETUP=""
GRIDOPTS=""
BKGVERR=""
ANBKGERR=""
JCOPTS=""
STRONGOPTS=""
OBSQC=""
OBSINPUT=""
SUPERRAD=""
LAGDATA=""
HYBRIDENSEMBLE=""
RR_CLDSURF=""
CHEM=""
SINGLEOB=""
SETUP_NSST=""
NSST=""


# T254 quadratic grid dimensions are nlon=768, nlat=384
# T254 linear grid dimensions are nlon=512, nlat=256+2
# HYBRIDENSEMBLE namelist below set up for linear T254 grid
if [[ "$DOHYBVAR" = "NO" ]]; then
   STRONGOPTS="tlnmc_option=1,"
fi
if [[ "$DOHYBVAR" = "YES" ]]; then
ensemble_dir='./ensemble_data/'
export ensemble_path=${ensemble_dir:-./}
HYBRIDENSEMBLE="l_hyb_ens=.true.,n_ens=$ENS_NUM_ANAL,beta_s0=0.25,readin_beta=.false.,s_ens_h=800.,s_ens_v=-0.8,generate_ens=.false.,uv_hyb_ens=.true.,jcap_ens=574,nlat_ens=578,nlon_ens=1152,aniso_a_en=.false.,jcap_ens_test=574,readin_localization=.true.,oz_univ_static=.false.,ensemble_path='${ensemble_path}',"
## ens_fast_read=.true.,"
fi

# Turn off generation of diagnostic files for GFS run
SETUPGFS=""
if [[ "$dumpobs" = "gfs" ]]; then
  SETUPGFS="diag_rad=.false.,diag_pcp=.false.,diag_conv=.false.,diag_ozone=.false.,write_diag(3)=.false.,"
fi

##SETUP_NSST="tzr_qc=1,"
##NSST="nst_gsi=3,nstinfo=4,zsea1=0,zsea2=5,fac_dtl=1,fac_tsl=1,"

SETUP_4DVAR=""
JCOPTS_4DVAR=""
STRONGOPTS_4DVAR=""
if [[ "$DO4DENSVAR" = "YES" ]]; then
## SETUP_4DVAR="l4densvar=.true.,ens_nstarthr=3,nhr_obsbin=1,nhr_assimilation=6,lwrite4danl=.true."   # 1 hourly, 7 analysis
   SETUP_4DVAR="l4densvar=.true.,ens_nstarthr=3,nhr_obsbin=1,nhr_assimilation=6,lwrite4danl=.false.,"  # 1 hourly, 1 analysis
   JCOPTS_4DVAR="ljc4tlevs=.true.,"
   STRONGOPTS_4DVAR="tlnmc_option=3,"
##   if [[ "$dumpobs" = "gfs" ]]; then
##    SETUP_4DVAR="l4densvar=.true.,ens4d_nstarthr=3,nhr_obsbin=1,nhr_assimilation=6,"
##   fi
fi

##SETUPLIMQ="factqmax=0.0,"
SETUPLIMQ=""

SETUP="$SETUP_4DVAR $SETUP_NSST $SETUPGFS $SETUPLIMQ"
JCOPTS="$JCOPTS_4DVAR"
STRONGOPTS="$STRONGOPTS_4DVAR"

export lrun_subdirs=${lrun_subdirs:-".true."}
##export use_readin_anl_sfcmask=".true."
export use_readin_anl_sfcmask=${use_readin_anl_sfcmask:-".false."}
export use_gfs_nemsio=".true."
export use_gfs_nemsio=${use_gfs_nemsio:-".false."}
crtm_coeffs=./crtm_coeffs/
export crtm_coeffs=${crtm_coeffs:-./}

# MLS ozone version 2.0 before 2013010306.   Version 3.0 starting 2013010306.
# Overwrite MLS entry in OBSINPUT for dates before 2013010306.
# NOTE:  scripting below DOES NOT WORK for GSI
##if [ $CDATE -lt 2013010306 ] ; then
##  export OBS_INPUT="mlsbufr        mls20       aura        mls20_aura          0.0     0     0"
##elif
##  export OBS_INPUT="mlsbufr        mls30       aura        mls30_aura          0.0     0     0"
##fi

rm -f gsiparm.anl

## use_readin_anl_sfcmask=${use_readin_anl_sfcmask},
cat <<EOF > gsiparm.anl
 &SETUP
   miter=1,niter(1)=5,niter(2)=150,
   niter_no_qc(1)=25,niter_no_qc(2)=0,
   write_diag(1)=.true.,write_diag(2)=.false.,write_diag(3)=.true.,
   qoption=2,
   gencode=82,factqmin=5.0,factqmax=0.005,deltim=$DELTIM,
   iguess=-1,
   oneobtest=.false.,retrieval=.false.,l_foto=.false.,
   use_pbl=.false.,use_compress=.true.,nsig_ext=12,gpstop=50.,
   use_gfs_nemsio=${use_gfs_nemsio},lrun_subdirs=${lrun_subdirs},use_readin_anl_sfcmask=${use_readin_anl_sfcmask},
   crtm_coeffs_path='${crtm_coeffs}',
   newpc4pred=.true.,adp_anglebc=.true.,angord=4,passive_bc=.true.,use_edges=.false.,
   diag_precon=.true.,step_start=1.e-3,emiss_bc=.true.,thin4d=.true.,cwoption=3,
   $SETUP
 /
 &GRIDOPTS
   JCAP_B=$JCAP_B,JCAP=$JCAP,NLAT=$NLAT,NLON=$NLON,nsig=$LEVS,
   regional=.false.,nlayers(63)=3,nlayers(64)=6,
   $GRIDOPTS
 /
 &BKGERR
   vs=0.7,
   hzscl=1.7,0.8,0.5,
   hswgt=0.45,0.3,0.25,
   bw=0.0,norsp=4,
   bkgv_flowdep=.true.,bkgv_rewgtfct=1.5,
   bkgv_write=.false.,
   cwcoveqqcov=.false.,
   $BKGERR   
 /
 &ANBKGERR
   anisotropic=.false.,
   $ANBKGERR
 /
 &JCOPTS
   ljcdfi=.false.,alphajc=0.0,ljcpdry=.true.,bamp_jcpdry=5.0e7,
   $JCOPTS   
 /
 &STRONGOPTS
   tlnmc_option=2,nstrong=1,nvmodes_keep=8,period_max=6.,period_width=1.5,
   baldiag_full=.false.,baldiag_inc=.false.,
   $STRONGOPTS   
 /
 &OBSQC
   dfact=0.75,dfact1=3.0,noiqc=.true.,oberrflg=.false.,c_varqc=0.04,
   use_poq7=.true.,qc_noirjaco3_pole=.true.,vqc=.true.,
   aircraft_t_bc=.true.,biaspredt=1000.0,upd_aircraft=.true.,
   $OBSQC
 /
 &OBS_INPUT
   dmesh(1)=145.0,dmesh(2)=150.0,dmesh(3)=100.0,time_window_max=3.0,
   $OBSINPUT
 /
OBS_INPUT::
!  dfile          dtype       dplat       dsis                dval    dthin dsfcalc
   prepbufr       ps          null        ps                  0.0     0     0
   prepbufr       t           null        t                   0.0     0     0
   prepbufr_profl t           null        t                   0.0     0     0
   prepbufr       q           null        q                   0.0     0     0
   prepbufr_profl q           null        q                   0.0     0     0
   prepbufr       pw          null        pw                  0.0     0     0
   prepbufr       uv          null        uv                  0.0     0     0
   prepbufr_profl uv          null        uv                  0.0     0     0
   satwndbufr     uv          null        uv                  0.0     0     0
   prepbufr       spd         null        spd                 0.0     0     0
   prepbufr       dw          null        dw                  0.0     0     0
   radarbufr      rw          null        rw                  0.0     0     0
   nsstbufr       sst         nsst        sst                 0.0     0     0
   gpsrobufr      gps_bnd     null        gps                 0.0     0     0
   ssmirrbufr     pcp_ssmi    dmsp        pcp_ssmi            0.0    -1     0
   tmirrbufr      pcp_tmi     trmm        pcp_tmi             0.0    -1     0
   sbuvbufr       sbuv2       n16         sbuv8_n16           0.0     0     0
   sbuvbufr       sbuv2       n17         sbuv8_n17           0.0     0     0
   sbuvbufr       sbuv2       n18         sbuv8_n18           0.0     0     0
   hirs3bufr      hirs3       n17         hirs3_n17           0.0     1     0
   hirs4bufr      hirs4       metop-a     hirs4_metop-a       0.0     1     1
   gimgrbufr      goes_img    g11         imgr_g11            0.0     1     0
   gimgrbufr      goes_img    g12         imgr_g12            0.0     1     0
   airsbufr       airs        aqua        airs_aqua           0.0     1     1
   amsuabufr      amsua       n15         amsua_n15           0.0     1     1
   amsuabufr      amsua       n18         amsua_n18           0.0     1     1
   amsuabufr      amsua       metop-a     amsua_metop-a       0.0     1     1
   airsbufr       amsua       aqua        amsua_aqua          0.0     1     1
   amsubbufr      amsub       n17         amsub_n17           0.0     1     1
   mhsbufr        mhs         n18         mhs_n18             0.0     1     1
   mhsbufr        mhs         metop-a     mhs_metop-a         0.0     1     1
   ssmitbufr      ssmi        f15         ssmi_f15            0.0     1     0
   amsrebufr      amsre_low   aqua        amsre_aqua          0.0     1     0
   amsrebufr      amsre_mid   aqua        amsre_aqua          0.0     1     0
   amsrebufr      amsre_hig   aqua        amsre_aqua          0.0     1     0
   ssmisbufr      ssmis       f16         ssmis_f16           0.0     1     0
   ssmisbufr      ssmis       f17         ssmis_f17           0.0     1     0
   ssmisbufr      ssmis       f18         ssmis_f18           0.0     1     0
   ssmisbufr      ssmis       f19         ssmis_f19           0.0     1     0
   gsnd1bufr      sndrd1      g12         sndrD1_g12          0.0     1     0
   gsnd1bufr      sndrd2      g12         sndrD2_g12          0.0     1     0
   gsnd1bufr      sndrd3      g12         sndrD3_g12          0.0     1     0
   gsnd1bufr      sndrd4      g12         sndrD4_g12          0.0     1     0
   gsnd1bufr      sndrd1      g11         sndrD1_g11          0.0     1     0
   gsnd1bufr      sndrd2      g11         sndrD2_g11          0.0     1     0
   gsnd1bufr      sndrd3      g11         sndrD3_g11          0.0     1     0
   gsnd1bufr      sndrd4      g11         sndrD4_g11          0.0     1     0
   gsnd1bufr      sndrd1      g13         sndrD1_g13          0.0     1     0
   gsnd1bufr      sndrd2      g13         sndrD2_g13          0.0     1     0
   gsnd1bufr      sndrd3      g13         sndrD3_g13          0.0     1     0
   gsnd1bufr      sndrd4      g13         sndrD4_g13          0.0     1     0
   iasibufr       iasi        metop-a     iasi_metop-a        0.0     1     1
   gomebufr       gome        metop-a     gome_metop-a        0.0     2     0
   omibufr        omi         aura        omi_aura            0.0     2     0
   sbuvbufr       sbuv2       n19         sbuv8_n19           0.0     0     0
   hirs4bufr      hirs4       n19         hirs4_n19           0.0     1     1
   amsuabufr      amsua       n19         amsua_n19           0.0     1     1
   mhsbufr        mhs         n19         mhs_n19             0.0     1     1
   tcvitl         tcp         null        tcp                 0.0     0     0
   seviribufr     seviri      m08         seviri_m08          0.0     1     0
   seviribufr     seviri      m09         seviri_m09          0.0     1     0
   seviribufr     seviri      m10         seviri_m10          0.0     1     0
   hirs4bufr      hirs4       metop-b     hirs4_metop-b       0.0     1     1
   amsuabufr      amsua       metop-b     amsua_metop-b       0.0     1     1
   mhsbufr        mhs         metop-b     mhs_metop-b         0.0     1     1
   iasibufr       iasi        metop-b     iasi_metop-b        0.0     1     1
   gomebufr       gome        metop-b     gome_metop-b        0.0     2     0
   atmsbufr       atms        npp         atms_npp            0.0     1     0
   crisbufr       cris        npp         cris_npp            0.0     1     0
   crisfsbufr     cris-fsr    npp         cris-fsr_npp        0.0     1     0
   gsnd1bufr      sndrd1      g14         sndrD1_g14          0.0     1     0
   gsnd1bufr      sndrd2      g14         sndrD2_g14          0.0     1     0
   gsnd1bufr      sndrd3      g14         sndrD3_g14          0.0     1     0
   gsnd1bufr      sndrd4      g14         sndrD4_g14          0.0     1     0
   gsnd1bufr      sndrd1      g15         sndrD1_g15          0.0     1     0
   gsnd1bufr      sndrd2      g15         sndrD2_g15          0.0     1     0
   gsnd1bufr      sndrd3      g15         sndrD3_g15          0.0     1     0
   gsnd1bufr      sndrd4      g15         sndrD4_g15          0.0     1     0
   oscatbufr      uv          null        uv                  0.0     0     0
   mlsbufr        mls30       aura        mls30_aura          0.0     0     0
   avhambufr      avhrr       metop-a     avhrr3_metop-a      0.0     1     0
   avhpmbufr      avhrr       n18         avhrr3_n18          0.0     1     0
   amsr2bufr      amsr2       gcom-w1     amsr2_gcom-w1       0.0     3     0
   gmibufr        gmi         gpm         gmi_gpm             0.0     3     0
   saphirbufr     saphir      meghat      saphir_meghat       0.0     3     0
   ahibufr        ahi         himawari8   ahi_himawari8       0.0     3     0
   rapidscatbufr  uv          null        uv                  0.0     0     0
::
  &SUPEROB_RADAR
   $SUPERRAD
 /
  &LAG_DATA
   $LAGDATA
 /
  &HYBRID_ENSEMBLE
   $HYBRIDENSEMBLE
 /
  &RAPIDREFRESH_CLDSURF
   dfi_radar_latent_heat_time_period=30.0,
   $RR_CLDSURF
 /
  &CHEM
   $CHEM
 /
  &NST
   $NSST
 /
 &SINGLEOB_TEST
   maginnov=0.1,magoberr=0.1,oneob_type='t',
   oblat=45.,oblon=180.,obpres=1000.,obdattim=${adate},
   obhourset=0.,
   $SINGLEOB
 /
EOF


# Set fixed files
#   berror   = forecast model background error statistics
#   specoef  = CRTM spectral coefficients
#   trncoef  = CRTM transmittance coefficients
#   emiscoef = CRTM coefficients for IR sea surface emissivity model
#   aerocoef = CRTM coefficients for aerosol effects
#   cldcoef  = CRTM coefficients for cloud effects
#   satinfo  = text file with information about assimilation of brightness temperatures
#   satangl  = angle dependent bias correction file (fixed in time)
#   pcpinfo  = text file with information about assimilation of prepcipitation rates
#   ozinfo   = text file with information about assimilation of ozone data
#   errtable = text file with obs error for conventional data (optional)
#   convinfo = text file with information about assimilation of conventional data
#   bufrtable= text file ONLY needed for single obs test (oneobstest=.true.)
#   bftab_sst= bufr table for sst ONLY needed for sst retrieval (retrieval=.true.)
#   aeroinfo = text file with information about assimilation of aerosol data

anavinfo=$fixgsi/global_anavinfo.l${LEVS}.txt
berror=$fixgsi/Big_Endian/global_berror.l${LEVS}y${NLAT}.f77
locinfo=$fixgsi/global_hybens_info.l${LEVS}.txt
satinfo=$fixgsi/global_satinfo.txt
cldradinfo=$fixgsi/cloudy_radiance_info.txt
scaninfo=$fixgsi/global_scaninfo.txt
satangl=$fixgsi/global_satangbias.txt
pcpinfo=$fixgsi/global_pcpinfo.txt
ozinfo=$fixgsi/global_ozinfo.txt
convinfo=$fixgsi/global_convinfo.txt
insituinfo=$fixgsi/global_insituinfo.txt
errtable=$fixgsi/prepobs_errtable.global
aeroinfo=$fixgsi/global_aeroinfo.txt
atmsbeaminfo=$fixgsi/atms_beamwidth.txt

emiscoef_IRwater=$fixcrtm/Nalli.IRwater.EmisCoeff.bin
emiscoef_IRice=$fixcrtm/NPOESS.IRice.EmisCoeff.bin
emiscoef_IRland=$fixcrtm/NPOESS.IRland.EmisCoeff.bin
emiscoef_IRsnow=$fixcrtm/NPOESS.IRsnow.EmisCoeff.bin
emiscoef_VISice=$fixcrtm/NPOESS.VISice.EmisCoeff.bin
emiscoef_VISland=$fixcrtm/NPOESS.VISland.EmisCoeff.bin
emiscoef_VISsnow=$fixcrtm/NPOESS.VISsnow.EmisCoeff.bin
emiscoef_VISwater=$fixcrtm/NPOESS.VISwater.EmisCoeff.bin
emiscoef_MWwater=$fixcrtm/FASTEM6.MWwater.EmisCoeff.bin
aercoef=$fixcrtm/AerosolCoeff.bin
cldcoef=$fixcrtm/CloudCoeff.bin

# csem data source path
csem=/scratch4/NCEPDEV/jcsda/save/Ming.Chen/CRTM_CSEM/csem_v1.0.0
csemalg=$csem/csem_model.registor
csemfix=$csem/fix

# Only need this file for single obs test
bufrtable=$fixgsi/prepobs_prep.bufrtable

# Only need this file for sst retrieval
bftab_sst=$fixgsi/bufrtab.012

# Copy executable and fixed files to $DATA
$ncp1 $gsiexec      ./gsi.x

$ncp1 $anavinfo     ./anavinfo
$ncp1 $berror       ./berror_stats
$ncp1 $locinfo      ./hybens_info
$ncp1 $satinfo      ./satinfo
$ncp1 $cldradinfo   ./cloudy_radiance_info.txt
$ncp1 $scaninfo     ./scaninfo
$ncp1 $pcpinfo      ./pcpinfo
$ncp1 $ozinfo       ./ozinfo
$ncp1 $convinfo     ./convinfo
$ncp1 $insituinfo   ./insituinfo
$ncp1 $errtable     ./errtable
$ncp1 $aeroinfo     ./aeroinfo
$ncp1 $atmsbeaminfo ./atms_beamwidth.txt

$ncp1 $bufrtable    ./prepobs_prep.bufrtable
$ncp1 $bftab_sst    ./bftab_sstphr

# Copy CRTM coefficient files based on entries in satinfo file
mkdir -p ${crtm_coeffs}
for file in `awk '{if($1!~"!"){print $1}}' satinfo | sort | uniq` ;do
   $ncp1 $fixcrtm/${file}.SpcCoeff.bin ${crtm_coeffs}
   $ncp1 $fixcrtm/${file}.TauCoeff.bin ${crtm_coeffs}
done
$ncp1 $emiscoef_IRwater  ${crtm_coeffs}Nalli.IRwater.EmisCoeff.bin
$ncp1 $emiscoef_IRice    ${crtm_coeffs}NPOESS.IRice.EmisCoeff.bin
$ncp1 $emiscoef_IRsnow   ${crtm_coeffs}NPOESS.IRsnow.EmisCoeff.bin
$ncp1 $emiscoef_IRland   ${crtm_coeffs}NPOESS.IRland.EmisCoeff.bin
$ncp1 $emiscoef_VISice   ${crtm_coeffs}NPOESS.VISice.EmisCoeff.bin
$ncp1 $emiscoef_VISland  ${crtm_coeffs}NPOESS.VISland.EmisCoeff.bin
$ncp1 $emiscoef_VISsnow  ${crtm_coeffs}NPOESS.VISsnow.EmisCoeff.bin
$ncp1 $emiscoef_VISwater ${crtm_coeffs}NPOESS.VISwater.EmisCoeff.bin
$ncp1 $emiscoef_MWwater  ${crtm_coeffs}FASTEM6.MWwater.EmisCoeff.bin
$ncp1 $aercoef           ${crtm_coeffs}AerosolCoeff.bin
$ncp1 $cldcoef           ${crtm_coeffs}CloudCoeff.bin

#copy csem data
$ncp1 $csemalg          ./csem_model.registor
cp -fR   $csemfix        csem_coeffs

# Copy observational data to $DATA
$ncp $datobs/${prefix_obs}.prepbufr                ./prepbufr
$ncp $datobs/${prefix_obs}.prepbufr.acft_profiles  ./prepbufr_profl
$ncp $datobs/${prefix_obs}.gpsro.${suffix}         ./gpsrobufr
$ncp $datobs/${prefix_obs}.satwnd.${suffix}        ./satwndbufr
$ncp $datobs/${prefix_obs}.spssmi.${suffix}        ./ssmirrbufr
$ncp $datobs/${prefix_obs}.sptrmm.${suffix}        ./tmirrbufr
$ncp $datobs/${prefix_obs}.osbuv8.${suffix}        ./sbuvbufr
$ncp $datobs/${prefix_obs}.goesfv.${suffix}        ./gsnd1bufr
$ncp $datobs/${prefix_obs}.1bamua.${suffix}        ./amsuabufr
$ncp $datobs/${prefix_obs}.1bamub.${suffix}        ./amsubbufr
$ncp $datobs/${prefix_obs}.1bhrs2.${suffix}        ./hirs2bufr
$ncp $datobs/${prefix_obs}.1bhrs3.${suffix}        ./hirs3bufr
$ncp $datobs/${prefix_obs}.1bhrs4.${suffix}        ./hirs4bufr
$ncp $datobs/${prefix_obs}.1bmhs.${suffix}         ./mhsbufr
$ncp $datobs/${prefix_obs}.1bmsu.${suffix}         ./msubufr
$ncp $datobs/${prefix_obs}.airsev.${suffix}        ./airsbufr
$ncp $datobs/${prefix_obs}.sevcsr.${suffix}        ./seviribufr
$ncp $datobs/${prefix_obs}.mtiasi.${suffix}        ./iasibufr
$ncp $datobs/${prefix_obs}.ssmit.${suffix}         ./ssmitbufr
## $ncp $datobs/${prefix_obs}.amsre.${suffix}         ./amsrebufr
$ncp $datobs/${prefix_obs}.ssmisu.${suffix}        ./ssmisbufr
$ncp $datobs/${prefix_obs}.gome.${suffix}          ./gomebufr
$ncp $datobs/${prefix_obs}.omi.${suffix}           ./omibufr
$ncp $datobs/${prefix_obs}.mls.${suffix}           ./mlsbufr
$ncp $datobs/${prefix_obs}.eshrs3.${suffix}        ./hirs3bufrears
$ncp $datobs/${prefix_obs}.esamua.${suffix}        ./amsuabufrears
$ncp $datobs/${prefix_obs}.esamub.${suffix}        ./amsubbufrears
$ncp $datobs/${prefix_obs}.atms.${suffix}          ./atmsbufr
$ncp $datobs/${prefix_obs}.cris.${suffix}          ./crisbufr
$ncp $datobs/${prefix_obs}.syndata.tcvitals.tm00   ./tcvitl
## $ncp $datobs/${prefix_obs}.gmi1cr.${suffix}        ./gmibufr
$ncp $datobs/${prefix_obs}.saphir.${suffix}        ./saphirbufr
$ncp $datobs/${prefix_obs}.iasidb.${suffix}        ./iasibufr_db
$ncp $datobs/${prefix_obs}.crisdb.${suffix}        ./crisbufr_db
$ncp $datobs/${prefix_obs}.atmsdb.${suffix}        ./atmsbufr_db
$ncp $datobs/${prefix_obs}.hrs3db.${suffix}        ./hirs3bufr_db
$ncp $datobs/${prefix_obs}.amuadb.${suffix}        ./amsuabufr_db
$ncp $datobs/${prefix_obs}.amubdb.${suffix}        ./amsubbufr_db
$ncp $datobs/${prefix_obs}.esiasi.${suffix}        ./iasibufrears
$ncp $datobs/${prefix_obs}.eshrs3.${suffix}        ./hirs3bufrears
$ncp $datobs/${prefix_obs}.esamua.${suffix}        ./amsuabufrears
$ncp $datobs/${prefix_obs}.esamub.${suffix}        ./amsubbufrears
$ncp $datobs/${prefix_obs}.nsstbufr                ./nsstbufr
$ncp $datobs/${prefix_obs}.avcsam.${suffix}        ./avhambufr
$ncp $datobs/${prefix_obs}.avcspm.${suffix}        ./avhpmbufr

$ncp $datges/${prefix_tbc}.abias                   ./satbias_in
$ncp $datges/${prefix_tbc}.abias_pc                ./satbias_pc
$ncp $datges/${prefix_tbc}.abias_air               ./aircftbias_in

$ncp $datges/${prefix_sfc}.sfcf003.nemsio          ./sfcf03
$ncp $datges/${prefix_sfc}.sfcf006.nemsio          ./sfcf06
$ncp $datges/${prefix_sfc}.sfcf009.nemsio          ./sfcf09

$ncp $datges/${prefix_sfc}.nstf003.nemsio          ./nstf03
$ncp $datges/${prefix_sfc}.nstf006.nemsio          ./nstf06
$ncp $datges/${prefix_sfc}.nstf009.nemsio          ./nstf09

$ncp $datatm/${prefix_atm}.atmgm3.nemsio           ./sigf03
$ncp $datatm/${prefix_atm}.atmges.nemsio           ./sigf06
$ncp $datatm/${prefix_atm}.atmgp3.nemsio           ./sigf09

$ncp $datobs/${prefix_obs}.sfcgcy                  ./sfcgcy

if [[ "$DO4DENSVAR" = "YES" ]]; then
   $ncp $datges/${prefix_sfc}.sfcf004.nemsio          ./sfcf04
   $ncp $datges/${prefix_sfc}.sfcf005.nemsio          ./sfcf05
   $ncp $datges/${prefix_sfc}.sfcf007.nemsio          ./sfcf07
   $ncp $datges/${prefix_sfc}.sfcf008.nemsio          ./sfcf08

   $ncp $datges/${prefix_sfc}.nstf004.nemsio          ./nstf04
   $ncp $datges/${prefix_sfc}.nstf005.nemsio          ./nstf05
   $ncp $datges/${prefix_sfc}.nstf007.nemsio          ./nstf07
   $ncp $datges/${prefix_sfc}.nstf008.nemsio          ./nstf08

   $ncp $datatm/${prefix_atm}.atmgm2.nemsio           ./sigf04
   $ncp $datatm/${prefix_atm}.atmgm1.nemsio           ./sigf05
   $ncp $datatm/${prefix_atm}.atmgp1.nemsio           ./sigf07
   $ncp $datatm/${prefix_atm}.atmgp2.nemsio           ./sigf08
fi

if [[ "$DOHYBVAR" = "YES" ]]; then
  flist="06"
  if [[ "$DO4DENSVAR" = "YES" ]]; then
     flist="03 04 05 06 07 08 09"
  fi
  mkdir -p $ensemble_path
  for fh in $flist; do
     SIGGESENS=${prefix_ens}.atmf0${fh}s
     imem=1
     while [[ $imem -le $ENS_NUM_ANAL ]]; do
        member="mem"`printf %03i $imem`
        sigens=${SIGGESENS}.${member}.nemsio
        $ncp $datens/$sigens ${ensemble_path}sigf${fh}_ens_${member}
        (( imem = $imem + 1 ))
     done
  done   
fi


if [[ "${use_readin_anl_sfcmask}" = ".true." ]]; then
   $ncp $datens/${prefix_ens}.sfcf006.ensmean.nemsio ./sfcf06_anlgrid
fi


# If requested, copy and de-tar guess radstat file
$ncp $datges/${prefix_sfc}.radstat ./radstat.gdas

# If requested, copy and de-tar guess radstat file
if [[ $USE_CFP = YES ]]; then
   rm $DATA/unzip.sh
   rm $DATA/mp_unzip.sh
   set +x
cat <<\EOFunzip > unzip.sh
#!/bin/ksh
{ echo
 set -aux
 diag_file=$1
 fname=`echo $diag_file | cut -d'.' -f1`
 date=`echo $diag_file | cut -d'.' -f2`
 $UNCOMPRESS $diag_file
 fnameges=$(echo $fname|sed 's/_ges//g')
 mv $fname.$date $fnameges
}
EOFunzip
   set -x
   chmod 755 $DATA/unzip.sh
fi

listdiag=`tar xvf radstat.gdas | cut -d' ' -f2 | grep _ges`
for type in $listdiag; do
   diag_file=`echo $type | cut -d',' -f1`
   if [[ $USE_CFP = YES ]] ; then
      echo "$DATA/unzip.sh $diag_file" | tee -a $DATA/mp_unzip.sh
   else
      fname=`echo $diag_file | cut -d'.' -f1`
      date=`echo $diag_file | cut -d'.' -f2`
      $UNCOMPRESS $diag_file
      fnameges=$(echo $fname|sed 's/_ges//g')
      mv $fname.$date $fnameges
   fi
done

if [[ $USE_CFP = YES ]] ; then
   chmod 755 $DATA/mp_unzip.sh
   ncmd=`cat $DATA/mp_unzip.sh | $wc -l`
   npe_node_a=24
   if [ $ncmd -lt $npe_node_a ]; then
      npe_node_a = $ncmd
   fi
   if [ $ncmd -gt 0 ]; then
      export APRUNCFP='aprun -q -j1 -n$ncmd -N$npe_node_a -d1 cfp'
      export APRUNCFP_UNZIP=$(eval echo $APRUNCFP)
      $APRUNCFP_UNZIP $DATA/mp_unzip.sh
   fi
fi



# Run gsi under Parallel Operating Environment (poe) on NCEP IBM

printenv > stdout.env
printenv

mpirun -np $PBS_NP $DATA/gsi.x < gsiparm.anl 1> stdout 2> stderr
rc=$?

if [[ "$rc" != "0" ]]; then
   echo "***ERROR*** gsi.x failed to run to completion, with an error code of $rc"
   exit
fi

if [[ "$GENDIAG" = "NO" ]] ; then
  exit
fi

# Save output
mkdir -p $SAVDIR
cat stdout fort.2* > $SAVDIR/stdout.anl.$adate
cat fort.2*        > $SAVDIR/gsistat.$dumpobs.$adate
$ncp1 siganl          $SAVDIR/gfnanl.$dumpobs.$adate
$ncp1 satbias_out     $SAVDIR/biascr.$dumpobs.$adate
$ncp1 satbias_pc.out  $SAVDIR/biascr_pc.$dumpobs.$adate
$ncp1 satbias_out.int $SAVDIR/biascr.int.$dumpobs.$adate

#$ncp1 sfcanl.gsi      $SAVDIR/sfcanl_gsi.$dumpobs.$adate
##$ncp1 sfcf06          $SAVDIR/sfcf06.$dumpges.$gdate
##$ncp1 sigf06          $SAVDIR/sigf06.$dumpges.$gdate


CNVSTAT=$SAVDIR/cnvstat.gdas.$adate
PCPSTAT=$SAVDIR/pcpstat.gdas.$adate
OZNSTAT=$SAVDIR/oznstat.gdas.$adate
RADSTAT=$SAVDIR/radstat.gdas.$adate

rm -f $CNVSTAT
rm -f $PCPSTAT
rm -f $OZNSTAT
rm -f $RADSTAT


cd $DATA    # we should already be in $DATA, but extra cd to be sure.
rm -rf diag_*


# Set up lists and variables for various types of diagnostic files.
ntype=3

diagtype[0]="conv"
diagtype[1]="pcp_ssmi_dmsp pcp_tmi_trmm"
diagtype[2]="sbuv2_n16 sbuv2_n17 sbuv2_n18 sbuv2_n19 gome_metop-a gome_metop-b omi_aura mls30_aura"
diagtype[3]="hirs2_n14 msu_n14 sndr_g08 sndr_g11 sndr_g12 sndr_g13 sndr_g08_prep sndr_g11_prep sndr_g12_prep sndr_g13_prep sndrd1_g11 sndrd2_g11 sndrd3_g11 sndrd4_g11 sndrd1_g12 sndrd2_g12 sndrd3_g12 sndrd4_g12 sndrd1_g13 sndrd2_g13 sndrd3_g13 sndrd4_g13 sndrd1_g14 sndrd2_g14 sndrd3_g14 sndrd4_g14 sndrd1_g15 sndrd2_g15 sndrd3_g15 sndrd4_g15 hirs3_n15 hirs3_n16 hirs3_n17 amsua_n15 amsua_n16 amsua_n17 amsub_n15 amsub_n16 amsub_n17 hsb_aqua airs_aqua amsua_aqua imgr_g08 imgr_g11 imgr_g12 imgr_g14 imgr_g15 ssmi_f13 ssmi_f15 hirs4_n18 hirs4_metop-a amsua_n18 amsua_metop-a mhs_n18 mhs_metop-a amsre_low_aqua amsre_mid_aqua amsre_hig_aqua ssmis_f16 ssmis_f17 ssmis_f18 ssmis_f19 ssmis_f20 iasi_metop-a hirs4_n19 amsua_n19 mhs_n19 seviri_m08 seviri_m09 seviri_m10 cris_npp cris-fsr_npp atms_npp hirs4_metop-b amsua_metop-b mhs_metop-b iasi_metop-b avhrr_n18 avhrr_metop-a amsr2_gcom-w1 gmi_gpm saphir_meghat ahi_himawari8"

diaglist[0]=listcnv
diaglist[1]=listpcp
diaglist[2]=listozn
diaglist[3]=listrad

diagfile[0]=$CNVSTAT
diagfile[1]=$PCPSTAT
diagfile[2]=$OZNSTAT
diagfile[3]=$RADSTAT

numfile[0]=0
numfile[1]=0
numfile[2]=0
numfile[3]=0


# Set diagnostic file prefix based on lrun_subdirs variable
if [ $lrun_subdirs = ".true." ]; then
   prefix=" dir.*/"
else
   prefix="pe*"
fi


# Collect diagnostic files as a function of loop and type.
loops="01 03"
for loop in $loops; do
   date
   case $loop in
      01) string=ges;;
      03) string=anl;;
       *) string=$loop;;
   esac
   n=-1
   while [ $((n+=1)) -le $ntype ] ;do
      for type in `echo ${diagtype[n]}`; do
         date
         count=`ls ${prefix}${type}_${loop}* | $wc -l`
         if [ $count -gt 0 ]; then
            cat ${prefix}${type}_${loop}* > diag_${type}_${string}.${CDATE}${DIAG_SUFFIX}
            echo "diag_${type}_${string}.${CDATE}*" >> ${diaglist[n]}
            numfile[n]=`expr ${numfile[n]} + 1`
         fi
      done
   done
done
date

cd $DATA    # we should already be in $DATA, but extra cd to be sure.

# If requested, compress diagnostic files
if [[ $DIAG_COMPRESS = YES ]]; then
   for file in `ls diag_*${CDATE}${DIAG_SUFFIX}`; do
      date
      $COMPRESS $file
   done
fi
date

# If requested, create diagnostic file tarballs
if [[ $DIAG_TARBALL = YES ]]; then
   n=-1
   while [ $((n+=1)) -le $ntype ] ;do
      date
      TAROPTS="-uvf"
      if [ ! -s ${diagfile[n]} ]; then
         TAROPTS="-cvf"
      fi
      if [ ${numfile[n]} -gt 0 ]; then
         tar $TAROPTS ${diagfile[n]} `cat ${diaglist[n]}`
      fi
   done

#  Restrict CNVSTAT 
   chmod 750 $CNVSTAT
   chgrp rstprod $CNVSTAT
fi
date


# If requested, clean up $DATA
if [[ "$CLEAN" = "YES" ]];then
   if [[ $rc -eq 0 ]];then
      rm -rf $DATA
      cd $DATA
      cd ../
      rmdir $DATA
   fi
fi
date

exit
