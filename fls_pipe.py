"""This module will take flc.fits files and apply Jay Anderson's selfcal fortran scripts,
outputing a final fls.fits file, which has edited SCI and DQ arrays.  Work to add some
kind of change to the ERR arrays has yet to be done.
"""


##LOCAL IMPORTS##
import wfcref

##GLOBAL IMPORTS##
import glob,subprocess,os,drizzlepac,shutil,time,stat
from astropy.io import fits
import numpy as NP


# ###### Run GALIGN #####   Assumes Input files are present, and correct
def galign(working_direc,ACSfilter):
    r"""Runs the first half of Jay's SelfCal scripts, does alignment
    of all images in folder.

    This function will take all the images in the working direcotry - 'working_direc'
    and runs them through the first half of Jay's codes, hst2galign.  The important
    products from this scrpt that are needed for the next step are the .mat files. A
    folder named 01.hst2galign will be created in your working directory where hst2galign
    will be run.

    Parameters
    ----------
    working_direc : string
        Working directory that contains the set of files you
        would like to process though hst2galign
    ACSfilter : string
        ACS filter of data in working_directory, case-insensitive
        i.e.  "F814W" or "f775w"

    Raises
    ------
    OSError
        If the function was unable to create a directory in 'working_dirc'

    See Also
    --------
    selfcal,final_fls

    Notes
    -----
    This function is using the hst2galign executable stored here:
    /internal/data1/frontier/software/selfcal/hst2galign.e

    Examples
    --------
    >>> import fls_pipe
    >>> fls_pipe.galign('working_direc/','f814w')
    """

    #time test
    start_time=time.time()
    
    #initial variable setup
    galign_direc=working_direc + '01.hst2galign'
    galign_exe='/internal/data1/frontier/software/selfcal/hst2galign.e'

    #directory setup and change directory
    if not os.path.exists(galign_direc): 
        try:
            os.mkdir(galign_direc,stat.S_IRWXU | stat.S_IRWXG | stat.S_IROTH | stat.S_IXOTH)
        except OSError as e:
            print("Oops '{}'".format(e.strerror()))
    
    #save current directory to return after galign is finished
    starting_direc=os.getcwd()

    #run galign
    os.chdir(galign_direc)
    if os.getcwd() == galign_direc:
        inst=' "INST=ACSWFC"'
        expt=' "EXPT_0=1307.0"'
        gcx=' "FILEGCX=/grp/webpages/jayander/GCLIB/ACS_WFC_SM4/wfc_'+ACSfilter.upper()+'_gcx_SM4.fits"'
        gcy=' "FILEGCY=/grp/webpages/jayander/GCLIB/ACS_WFC_SM4/wfc_'+ACSfilter.upper()+'_gcy_SM4.fits"'
        flc=' ' + working_direc + 'j*flc.fits'
        exe_call=galign_exe+inst+expt+gcx+gcy+flc
        subprocess.call(exe_call,shell=True)
        print('finished succesfully!, changing permissions')
        
        #change permissions
        CHlist=glob.glob('*')
        for file in CHlist:
            os.chmod(file,stat.S_IRWXU | stat.S_IRWXG | stat.S_IROTH | stat.S_IXOTH)
        
        #change back to original starting directory
        os.chdir(starting_direc)
    else:
        print 'could not change directory to hst2galign directory'

    #time test
    end_time=time.time()
    print end_time-start_time,'seconds'

# ###### Run SELFCAL #####
def selfcal(working_direc,ACSfilter):
    r"""Runs the second half of Jay's SelfCal scripts, creates final delta_dark file.

    This function will copy over the .mat files produced from hst2galign, run the
    second half of Jay's code, hst2selfcal on the images in the working directory
    and runs them through the first half of Jay's codes, hst2galign.  The important
    products from this scrpt that are needed for the next step are the 
    dark_bar.?.fits files. A folder named 02.hst2selfcal will be created in your
    working directory where hst2selfcal will be run.

    Parameters
    ----------
    working_direc : string
        Working directory that contains the set of files you
        would like to process though hst2selfcal
    ACSfilter : string
        ACS filter of data in working_directory, case insensitive
        i.e.  "F814W" or "f775w"

    Raises
    ------
    OSError
        If the function was unable to create a directory in 'working_dirc'

    See Also
    --------
    galign,final_fls

    Notes
    -----
    This function is using the hst2selfcal executable stored here:
    /internal/data1/frontier/software/selfcal/hst2selfcal.e

    Examples
    --------
    >>> import fls_pipe
    >>> fls_pipe.selfcal('working_direc/','f814w')
    """
    
    #time test
    start_time=time.time()

    #initial variable setup
    selfcal_direc=working_direc + '02.hst2selfcal'
    selfcal_exe='/internal/data1/frontier/software/selfcal/hst2selfcal.e'

    #directory setup
    if not os.path.exists(selfcal_direc):
        try:
            os.mkdir(selfcal_direc,stat.S_IRWXU | stat.S_IRWXG | stat.S_IROTH | stat.S_IXOTH)
        except OSError as e:
            print("Oops '{}'".format(e.strerror()))
        
    filelistx=glob.glob(working_direc+'01.hst2galign/*mat')
    for filex in filelistx:
        shutil.copy2(filex,selfcal_direc)

    #save current directory to return after selfcal is finished
    starting_direc=os.getcwd()

    #run selfcal
    os.chdir(selfcal_direc)
    if os.getcwd() == selfcal_direc:
        expt=' "EXPT_0=1307.0"'
        gcx=' "FILEGCX=/grp/webpages/jayander/GCLIB/ACS_WFC_SM4/wfc_'+ACSfilter.upper()+'_gcx_SM4.fits"'
        gcy=' "FILEGCY=/grp/webpages/jayander/GCLIB/ACS_WFC_SM4/wfc_'+ACSfilter.upper()+'_gcy_SM4.fits"'
        flc=' ' + working_direc + 'j*flc.fits'
        exe_call=selfcal_exe+expt+gcx+gcy+flc
        subprocess.call(exe_call,shell=True)
        print('finished succesfully!')

        #change permissions
        CHlist=glob.glob('*')
        for file in CHlist:
            os.chmod(file,stat.S_IRWXU | stat.S_IRWXG | stat.S_IROTH | stat.S_IXOTH)

        #change to original directory
        os.chdir(starting_direc)
    else:
        print 'could not change directory'

    #time test
    end_time=time.time()
    print end_time-start_time,'seconds'

# ######## Re-flag DQ array in temp prep file########
def final_fls(working_direc):
    r"""Applies the selfcal delta_dark from Jay's SelfCal scripts to flc files.

    This function will copy over the delta_dark file produced from hst2selfcal, and
    subtract it from the flc.fits files to make fls.fits files. The script will
    also re-inititalize the DQ array of the images after the delta_dark subtraction
    by adding the delta_dark to the superdark file, re-calcuating the DQ array, and
    changing these flags in the corresponding fls.fits file. A folder named 03.finalfls
    will be created in your working directory where the processing will take place.

    Parameters
    ----------
    working_direc : string
        Working directory that contains the set of files you
        would like to process into fls.fits files

    Raises
    ------
    OSError
        If the function was unable to create/remove  a directory in 'working_dirc'

    See Also
    --------
    galign,selfcal

    Notes
    -----
    

    Examples
    --------
    >>> import fls_pipe
    >>> fls_pipe.final_fls('working_direc/')
    """
    
    #time test
    start_time=time.time()

    #initial variable setup
    exptime=1307.0
    fls_direc=working_direc + '03.finalfls'
    selfcal_dark=working_direc + '02.hst2selfcal/dark_bar.03.fits'

    #save current directory to return after fls is finished
    starting_direc=os.getcwd()

    #directory setup, wipe directory if it already exists
    if os.path.exists(fls_direc):
        try:
            shutil.rmtree(fls_direc)
        except OSError as e:
            print("Oops '{}'".format(e.strerror()))
        #What other errros might come up?

    os.mkdir(fls_direc,stat.S_IRWXU | stat.S_IRWXG | stat.S_IROTH | stat.S_IXOTH)

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

        #record and clear out flag 256
        darkpp_hdul=fits.open(darkpp)
        mask1_chip2=NP.where(darkpp_hdul[3].data & 256, True, False)
        mask1_chip1=NP.where(darkpp_hdul[6].data & 256, True, False)
        darkpp_hdul.close()

        drizzlepac.resetbits.reset_dq_bits(darkpp,256) 

        #Run deltadark subtraction, and DQ array correction
        dqnew_hdul=fits.open(darkpp)
        for files in filelists:
            #record flag 256 from sci image
            fls_hdul=fits.open(files)
            mask2_chip2=NP.where(fls_hdul[3].data & 256, True, False)
            mask2_chip1=NP.where(fls_hdul[6].data & 256, True, False)
            fls_hdul.close()

            #clearn out DQ flags from sci image
            drizzlepac.resetbits.reset_dq_bits(files,16+64+256)
            
            #add DQ flags from deltadark into sci DQ array
            fls_hdul=fits.open(files,mode='update')
            fls_expt=fls_hdul[0].header['EXPTIME']
    
            fls_hdul[1].data[:,:] -= fls_expt*delta_hdul[0].data[:2048]/exptime
            fls_hdul[4].data[:,:] -= fls_expt*delta_hdul[0].data[2048:]/exptime
    
            fls_hdul[3].data[:,:] += dqnew_hdul[3].data[:,:]
            fls_hdul[6].data[:,:] += dqnew_hdul[6].data[:,:]

            #add back in DQ flag 256 from sci and deltadark frame, using array masks
            finmask_chip2=mask1_chip2 | mask2_chip2
            finmask_chip1=mask1_chip1 | mask2_chip1
            Marray_chip2=NP.ma.array(NP.zeros((2048,4096)),mask=finmask_chip2,dtype=int)
            Marray_chip1=NP.ma.array(NP.zeros((2048,4096)),mask=finmask_chip1,dtype=int)
            DQ256_chip2=Marray_chip2.filled(fill_value=256)
            DQ256_chip1=Marray_chip1.filled(fill_value=256)

            fls_hdul[3].data[:,:] += DQ256_chip2[:,:]
            fls_hdul[6].data[:,:] += DQ256_chip1[:,:]
    
            fls_hdul.close()
    
        dqnew_hdul.close()
        delta_hdul.close()

        #change permissions
        CHlist=glob.glob('*')
        for file in CHlist:
            os.chmod(file,stat.S_IRWXU | stat.S_IRWXG | stat.S_IROTH | stat.S_IXOTH)

        print('finished succesfully!')
        os.chdir(starting_direc)

    #time test
    end_time=time.time()
    print end_time-start_time,'seconds'
