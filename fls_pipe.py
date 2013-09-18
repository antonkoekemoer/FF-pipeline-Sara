##GLOBAL VARIABLE##
working_direc='/internal/data1/frontier/Soft_Testing/'
exptime = 1307.0

##LOCAL IMPORTS##
import wfcref

##GLOBAL IMPORTS##
import glob,subprocess,os,drizzlepac,shutil
from astropy.io import fits
import numpy as NP


# ###### Run GALIGN #####   Assumes Input files are present, and correct
def galign(working_direc):
    galign_direc=working_direc + '01.hst2galign'
    galign_exe='/internal/data1/frontier/software/selfcal/hst2galign.e'

    #directory setup and change directory
    if not os.path.exists(galign_direc): 
        try:
            os.mkdir(galign_direc)
        except OSError as e:
            print("Oops '{}'".format(e.strerror()))

    #run galign
    os.chdir(galign_direc)
    if os.getcwd() == galign_direc:
        inst=' "INST=ACSWFC"'
        expt=' "EXPT_0=1307.0"'
        gcx=' "FILEGCX=/grp/webpages/jayander/GCLIB/ACS_WFC_SM4/wfc_F814W_gcx_SM4.fits"'
        gcy=' "FILEGCY=/grp/webpages/jayander/GCLIB/ACS_WFC_SM4/wfc_F814W_gcy_SM4.fits"'
        flc=' ' + working_direc + 'j*flc.fits'
        exe_call=galign_exe+inst+expt+gcx+gcy+flc
        subprocess.call(exe_call,shell=True)
    else:
        print 'could not change directory to hst2galign directory'

# ###### Run SELFCAL #####
def selfcal(working_direc):
    selfcal_direc=working_direc + '02.hst2selfcal'
    selfcal_exe='/internal/data1/frontier/software/selfcal/hst2selfcal.e'

    #directory setup
    if not os.path.exists(selfcal_direc):
        try:
            os.mkdir(selfcal_direc)
        except OSError as e:
            print("Oops '{}'".format(e.strerror()))
        
    filelistx=glob.glob(working_direc+'01.hst2galign/*xym*')
    for filex in filelistx:
        shutil.copy2(filex,selfcal_direc)

    #run selfcal
    os.chdir(selfcal_direc)
    if os.getcwd() == selfcal_direc:
        expt=' "EXPT_0=1307.0"'
        gcx=' "FILEGCX=/grp/webpages/jayander/GCLIB/ACS_WFC_SM4/wfc_F814W_gcx_SM4.fits"'
        gcy=' "FILEGCY=/grp/webpages/jayander/GCLIB/ACS_WFC_SM4/wfc_F814W_gcy_SM4.fits"'
        flc=' ' + working_direc + 'j*flc.fits'
        exe_call=selfcal_exe+expt+gcx+gcy+flc
        subprocess.call(exe_call,shell=True)
    else:
        print 'could not change directory'


# ######## Re-flag DQ array in temp prep file########
def final_fls(working_direc,exptime):
    fls_direc=working_direc + '03.finalfls'
    selfcal_dark=working_direc + '02.hst2selfcal/dark_bar.03.fits'

    #directory setup, wipe directory if it already exists
    if os.path.exists(fls_direc):
        try:
            shutil.rmtree(fls_direc)
        except OSError as e:
            print("Oops '{}'".format(e.strerror()))
        #What other errros might come up?

    os.mkdir(fls_direc)

    shutil.copy2(selfcal_dark,fls_direc)
    filelistc=glob.glob(working_direc+'*flc.fits')
    for filec in filelistc:
        new_fname=fls_direc + filec[-19:].replace('flc','fls')
        shutil.copy2(filec,new_fname)
    
    #run final processing
    os.chdir(fls_direc)
    if os.getcwd() == fls_direc:
        #copy over dummy dark
        filelists=glob.glob('*fls.fits') 
        darkpp=fits.getval(filelists[0],'DARKFILE')[5:]
        darkpp_loc='/grp/hst/cdbs/jref/'+darkpp
        shutil.copy2(darkpp_loc,'.')

        #zero out DQ array in dummy dark
        wfcref.zero_DQ_arrays(darkpp)

        #add selfcal to dummy dark and rerun DQ array flagging
        darkpp_hdul=fits.open(darkpp,mode='update')
        delta_hdul=fits.open('dark_bar.03.fits')
        darkpp_hdul[1].data[:,:] += delta_hdul[0].data[:2048,:]/exptime
        darkpp_hdul[4].data[:,:] += delta_hdul[0].data[2048:,:]/exptime
        darkpp_hdul.close()

        wfcref.hot_pixels(darkpp,'acs', do_warm=1, do_sat=1)

        #Run deltadark subtraction, and DQ array correction
        dqnew_hdul=fits.open(darkpp)
        for files in filelists:
            drizzlepac.resetbits.reset_dq_bits(files,16+64+256)
            fls_hdul=fits.open(files,mode='update')
            fls_expt=fls_hdul[0].header['EXPTIME']
    
            fls_hdul[1].data[:,:] -= fls_expt*delta_hdul[0].data[:2048]/exptime
            fls_hdul[4].data[:,:] -= fls_expt*delta_hdul[0].data[2048:]/exptime
    
            fls_hdul[3].data[:,:] += dqnew_hdul[3].data[:,:]
            fls_hdul[6].data[:,:] += dqnew_hdul[6].data[:,:]
    
            fls_hdul.close()
    
        dqnew_hdul.close()
        delta_hdul.close()
