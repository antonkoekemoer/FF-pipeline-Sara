"""This module will take flc.fits files and apply Jay Anderson's selfcal fortran scripts,
outputing a final fls.fits file, which has edited SCI and DQ arrays.  Work to add some
kind of change to the ERR arrays has yet to be done.
"""


##LOCAL IMPORTS##
import wfcref

##GLOBAL IMPORTS##
import glob,subprocess,os,drizzlepac,shutil
from astropy.io import fits


# ###### Run GALIGN #####   Assumes Input files are present, and correct
def galign(working_direc):
    r"""Runs the first half of Jay's SelfCal scripts, does alignment
    of all images in folder.

    This function will take all the images in the working direcotry - 'working_direc'
    and runs them through the first half of Jay's codes, hst2galign.  The important
    products from this scrpt that are needed for the next step are the .xym files. A
    folder named 01.hst2galign will be created in your working directory where hst2galign
    will be run.

    Parameters
    ----------
    working_direc : string
        Working directory that contains the set of files you
        would like to process though hst2galign

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
    >>> fls_pipe.galign('working_direc/')
    """
    
    #initial variable setup
    galign_direc=working_direc + '01.hst2galign'
    galign_exe='/internal/data1/frontier/software/selfcal/hst2galign.e'

    #directory setup and change directory
    if not os.path.exists(galign_direc): 
        try:
            os.mkdir(galign_direc)
        except OSError as e:
            print("Oops '{}'".format(e.strerror()))
    
    #save current directory to return after galign is finished
    starting_direc=os.getcwd()

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
        print('finished succesfully!')
        os.chdir(starting_direc)
    else:
        print 'could not change directory to hst2galign directory'

# ###### Run SELFCAL #####
def selfcal(working_direc):
    r"""Runs the second half of Jay's SelfCal scripts, creates final delta_dark file.

    This function will copy over the .xym files produced from hst2galign, run the
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
    >>> fls_pipe.selfcal('working_direc/')
    """
    
    #initial variable setup
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

    #save current directory to return after selfcal is finished
    starting_direc=os.getcwd()

    #run selfcal
    os.chdir(selfcal_direc)
    if os.getcwd() == selfcal_direc:
        expt=' "EXPT_0=1307.0"'
        gcx=' "FILEGCX=/grp/webpages/jayander/GCLIB/ACS_WFC_SM4/wfc_F814W_gcx_SM4.fits"'
        gcy=' "FILEGCY=/grp/webpages/jayander/GCLIB/ACS_WFC_SM4/wfc_F814W_gcy_SM4.fits"'
        flc=' ' + working_direc + 'j*flc.fits'
        exe_call=selfcal_exe+expt+gcx+gcy+flc
        subprocess.call(exe_call,shell=True)
        print('finished succesfully!')
        os.chdir(starting_direc)
    else:
        print 'could not change directory'


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

    print('finished succesfully!')
    os.chdir(starting_direc)
