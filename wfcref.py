"""
Collection of routines called by `wfc_reference.py`.

:Authors: Pey Lian Lim (ACS), A. Martel and T. Borders (WFC3)

:Organization: Space Telescope Science Institute

:History:
    * 2012/08/14 Ogaz trimmed down headers for ACS dark and bias frames
    * 2012/03/15 PLL added support for CALACS 8.0 that produces everything in ELECTRONS, not DN.
    * 2012/03/01 PLL changed CDRKFILE to DRKCFILE.
    * 2011/08/30 PLL modified `combine` and added `crx2out`.
    * 2011/07/29 PLL removed IRAF dependencies. Changed documentation format to Sphinx.
    * 2010/03/26 PLL changed `thres_bia` from 7.4 to 8.
    * 2010/03/04 PLL modified documentation to comply with Epydoc convention.
    * 2009/12/18 PLL fixed iref/jref selection bug in `electrons_normalize`.
    * 2009/12/17 PLL added anneal date for DESCRIP in `cdbs_update_header`.
    * 2009/12/14 PLL added subarray logic for ACS superbias. Improved algorithms.
    * 2009/12/09 PLL changed `thres_bia` from 8 to 7.4.
    * 2009/12/08 PLL added more header info to  comply with CDBS.
    * 2009/12/02 PLL added bias jump correction codes.
    * 2009/12/01 PLL changed some IRAF `imcalc` operations to PyFITS.
    * 2009/11/24 PLL added DQ flagging logic for ACS.
    * 2009/11/17 PLL obtained the codes written by AM from TMB. Conversion for ACS WFC started. Added `biasFile`, `saveTemp`, and `calname`.
    * 2-13/04/29 Ogaz changed pyfits to stpyfits
    
"""

# External modules
import os, subprocess, shutil, glob, time, numpy
from stsci.tools import stpyfits as pyfits
from stsci.tools import asnutil, parseinput
from stsci.convolve import boxcar

# Local modules
import mysky

#------- File type dependent functions. -------

def after_bias_combine(outputfinal, doSubArr=False):
    """ 
    The `calacs` output is the sum of the input RAW
    bias frames. So, divide by the number of input frames
    to get the final, mean reference bias frame. The
    integration time is 0 sec.

    Parameters
    ----------
    outputfinal: string
        The final, cr-rejected reference file.
        
    doSubArr: bool
        Use logic for subarray instead of full-frame.

    """
    print '  Calculating the mean bias frame (SCI and ERR).'

    if doSubArr:
        exts = [1, 2]
    else:
        exts = [1, 2, 4, 5]

    pf = pyfits.open(outputfinal, mode='update')
    c  = pf[0].header['CRSPLIT']
    for ext in exts: pf[ext].data /= c
    pf.close()

def bias_struct(outputfinal, instrument, saveTemp, doSubArr=False):
    """
    Identify bias structures and flag them as 128.
    ACS logic was obtained from M. Mutchler.

    .. note:: Flags are additive. WFC3 does not currently do this.

    Parameters
    ----------
    outputfinal: string
        Superbias image to flag.
    
    instrument: {'acs', 'wfc3'}.
        Instrument name.
    
    saveTemp: bool
        Keep temp smoothed bias images?

    doSubArr: bool
        Use logic for subarray instead of full-frame.

    """
    ext0, ext1 = 0, 1
    pf = pyfits.open(outputfinal, mode='update')

    # Properties for subarray
    if doSubArr:
        # Nothing is defined for WFC3 subarray
        if instrument == 'wfc3':
            print 'Flagging for DQ extension is skipped.'
            pf.close()
            return
        
        # Extensions
        ext_sci = [1]
        ext_dqa = [3]
        n_phys  = 22  # Num of pix in physical overscan

        # Quadrant coords excluding dummy overscans
        ccdamp = pf[ext0].header['CCDAMP']
        xx = pf[ext1].header['NAXIS1']
        yy = pf[ext1].header['NAXIS2']
        if ccdamp == 'A' or ccdamp == 'C':
            crd = [ [ [n_phys, xx, 0, yy] ] ]   # (x1 x2 y1 y2)
        elif ccdamp == 'B' or ccdamp == 'D':
            crd = [ [ [0, xx-n_phys, 0, yy] ] ] # (x1 x2 y1 y2)
        else:
            print 'Flagging for DQ extension is skipped for subarray AMP', ccdamp
            pf.close()
            return

        # Amps in quadrants
        amps = [ [ ccdamp ] ]
        regr = range(0, 1)  # 1 quadrant per ext

    # Properties for fullframe
    else:
        # Extensions
        ext_sci = [1, 4]
        ext_dqa = [3, 6]

        # Quadrant coords excluding dummy overscans
        if instrument == 'wfc3':
            crd = [ [ [25, 2073,  0, 2051], [2133, 4181,  0, 2051] ],   # Ext 1: left right (x1 x2 y1 y2)
                    [ [25, 2073, 19, 2070], [2133, 4181, 19, 2070] ] ]  # Ext 4: left right (x1 x2 y1 y2)
        else: # acs
            crd = [ [ [24, 2072,  0, 2048], [2072, 4120,  0, 2048] ],   # Ext 1: left right (x1 x2 y1 y2)
                    [ [24, 2072, 20, 2068], [2072, 4120, 20, 2068] ] ]  # Ext 4: left right (x1 x2 y1 y2)

        # Amps in quadrants
        amps = [ ['C', 'D'], ['A', 'B']]
        regr = range(0, 2)  # 2 quadrant per ext

    # Define smoothing box (smear more in X)
    box_x  = 11
    box_y  = 5
    box_sz = (box_y, box_x)

    # Define flag (see data handbook)
    flag_bia = 128
    thres_bia = 8.0 # Max increased to this value post-SM4.

    # Process each quadrant without dummy overscans   
    for i,sci in enumerate(ext_sci):
        dqa = ext_dqa[i]
        for j in regr:
            # Extract SCI for quadrant (y1 y2 x1 x2)
            im_sci = pf[sci].data[ crd[i][j][2]:crd[i][j][3], crd[i][j][0]:crd[i][j][1] ]
            im_sz  = im_sci.size

            # Smooth SCI data for flagging
            im_smo = boxcar(im_sci, box_sz)

            # Find pixels above threshold in smoothed image
            idx = numpy.where( im_sci >= (im_smo + thres_bia) )
            n   = len(idx[0])
            pct = n * 100.0 / im_sz
            print '... %i/%i pixels flagged in Amp %s (%.3f%%)' % (n, im_sz, amps[i][j], pct)
            if n < 1: continue

            # Flag DQ for quadrant
            im_dqa = pf[dqa].data[ crd[i][j][2]:crd[i][j][3], crd[i][j][0]:crd[i][j][1] ]
            im_dqa[idx] += flag_bia

    # Add HISTORY
    pf[ext0].header.add_history('Bias struct >= smoothed + %3.1f. Smooth window %ix%i.' % (thres_bia, box_x, box_y))
    pf.close()

def bia_for_quad(outputfinal, rootName, doSubArr=False):
    """
    Make a copy of superbias for WFC3 quad filters.

    .. note:: The output is just a copy with keyword changes to CCDAMP and APERTURE in order to accommodate revised CDBS selection rules. This output file will have a '_quad_bia.fits' attached to it. This is not used for ACS.

    Parameters
    ----------
    outputfinal: string
        Superbias to copy.

    rootName: string
        Prefix used to generate new file name.

    doSubArr: bool
        Use logic for subarray instead of full-frame.

    """

    # Do not proceed if input is not WFC3 full-frame
    if doSubArr:
        print 'Superbias copy for WFC3 quad filters is not needed for', outputfinal
        return

    ext0, ext1 = 0, 1

    outputfinal_cp = rootName.split('_')[0] + '_quad_bia.fits'
    shutil.copyfile(outputfinal, outputfinal_cp)
    
    fd = pyfits.open(outputfinal_cp, mode='update')
    
    fd[ext0].header.update('DESCRIP', pad_dscrp('Bias intended for use with quad filter and corner subarrays.'))
    fd[ext0].header.update('CCDAMP', 'QUAD_CORNER_SUBARRAYS')
    fd[ext0].header.update('APERTURE', 'SINGLE_AMP')

    fd.close()

def hot_pixels(outputfinal, instrument, do_warm=0, do_sat=0):
    """
    Identify the hot pixels in the science extensions
    (1 and 4) and flag them in the DQ extensions (3 and 6)
    with flag 16. Note that the flags are additive in the
    DQ arrays. The `outputfinal` image has unit of e-/sec.

    The hot pixels are defined above some lower threshold.

    .. note:: Flags are additive. WFC3 need to define new warm and saturation thresholds or never turn on `do_warm` and `do_sat` flags.

    **WFC3**

    For ground data, we pick 60 e-/hour in 1x1, -82 C,
    ABCD, gain=1.5 dark frames. This value is based on
    the tail of the histograms from a first pass of the
    superdarks. It's very conservative (could go as low
    as 40 e-/hour). This threshold is always the same
    for all binnings i.e. do NOT scale by factors of 4,
    9, 16 for binnings of 2, 3, and 4.
    
    For the warm ground darks (-49C), we picked a
    threshold of 4000 e-/hour (the average level is
    about 1500 e-/hour and the hot pixels are about
    >4500 e-/hour).

    **ACS**

    Thresholds are from M. Mutchler. The number for
    saturation was historical and meant to be used
    with 1000s DARK with gain 2. Hot and warm pix
    values are arbitrary.

    Parameters
    ----------
    outputfinal: string
        Superdark image to flag.
    
    instrument: {'acs', 'wfc3'}
        Instrument name.
    
    do_warm: {0, 1}
        Do warm pixel flagging?
            * 0: No
            * 1: Yes

    do_sat: {0, 1}
        Do saturated pixel flagging?
            * 0: No
            * 1: Yes

    """

    # Extensions for full-frame
    ext_sci = [1, 4]
    ext_dqa = [3, 6]
    ext0 = 0

    # Define flags (see data handbook)
    flag_hpx = 16
    flag_wpx = 64
    flag_sat = 256

    # Define flagging thresholds (e-/sec/pix)

    sec_in_hr = 3600.0

    if instrument == 'wfc3':

        threshold_hdr = 60.0 # e-/hr/pix
        threshold_unit = 'e-/hour'
        
        hot_thres_lo = threshold_hdr / sec_in_hr # hot pix
        
        # The rest of these are undefined for WFC3.
        # Using ACS values for now.
        wrm_thres_lo = 0.04   # warm pix
        sat_thres_lo = 60.0   # saturation
        
    else: # acs

        hot_thres_lo = 0.08   # hot pix
        wrm_thres_lo = 0.04   # warm pix
        sat_thres_lo = 60.0   # saturation

        threshold_hdr = hot_thres_lo
        threshold_unit = 'e-/sec'

    # Flag DQ arrays

    fd = pyfits.open(outputfinal, mode='update')

    for i, sci in enumerate(ext_sci):
        dqa = ext_dqa[i]
        im_sci = fd[sci].data
        im_dqa = fd[dqa].data
        im_fac = 100.0 / im_sci.size

        # Flag hot pixels
        idx_hot = numpy.where( im_sci >= hot_thres_lo)
        n_hot = len(idx_hot[0])
        p_hot = n_hot * im_fac
        print '... %i/%i hot pixels flagged as %i for EXT %i (%.3f%%)' % (n_hot, im_sci.size, flag_hpx, sci, p_hot)
        if n_hot < 1: continue
        im_dqa[idx_hot] += flag_hpx
        
        # Flag warm pixels
        if do_warm:
            idx_wrm = numpy.where( (im_sci >= wrm_thres_lo) & (im_sci < hot_thres_lo))
            n_wrm = len(idx_wrm[0])
            p_wrm = n_wrm * im_fac
            print '... %i/%i warm pixels flagged as %i for EXT %i (%.3f%%)' % (n_wrm, im_sci.size, flag_wpx, sci, p_wrm)
            if n_wrm < 1: continue
            im_dqa[idx_wrm] += flag_wpx

        # Flag saturated pixels
        if do_sat:
            idx_sat = numpy.where( im_sci >= sat_thres_lo)
            n_sat = len(idx_sat[0])
            p_sat = n_sat * im_fac
            print '... %i/%i saturated pixels flagged as %i for EXT %i (%.3f%%)' % (n_sat, im_sci.size, flag_sat, sci, p_sat)
            if n_sat < 1: continue
            im_dqa[idx_sat] += flag_sat

    # Add HISTORY

    fd[ext0].header.add_history('Hot pixel lower limit: %8.3f %s' % (threshold_hdr, threshold_unit))
    if do_warm: fd[ext0].header.add_history('Warm pix >= %5.3f and < %5.3f %s' % (wrm_thres_lo, hot_thres_lo, 'e-/sec'))
    if do_sat: fd[ext0].header.add_history('Saturated pix >= %5.1f %s' % (sat_thres_lo, 'e-/sec'))

    fd.close()

def amp_jump_corr(outputfinal, instrument, doplot=0):
    """
    Perform correction for bias jump across amp.

    Evaluate and remove the residual amplifier jump from
    base-dark frames. The estimate of the signal level
    in both side of the amplifier edge is done with
    fitting to the histogram, using the mode. Jump is
    removed only if above a certain threshold.

    .. note:: This was written by PLL based on `ampjump.cl` by M. Sirianni. Current threshold is 7.7e-6. PLL is not sure about this number. Currently not used by WFC3.

    Parameters
    ----------
    outputfinal: string
        Image to be corrected.
    
    instrument: {'acs', 'wfc3'}
        Currently not used in algorithm. For future.
    
    doplot: {0, 1}
        Display fits and verbose output? For testing.
            * 0: No
            * 1: Yes

    """

    #-------------------------
    # Define threshold.
    # Remove the jump only if it is > 0.01 e- in a 1300s exposure ~/half
    # orbit (7.7e-6 in 1 sec).

    # !!!!!!!!!!!!!!!!!!!!!!!
    # PLL: Not very sure about this number.
    #      It always trigger a correction with test data.
    # !!!!!!!!!!!!!!!!!!!!!!!
    
    threshold = 7.7e-6
        
    #-------------------------
    # Define extensions and quadrants.

    ext_sci = [1, 4]
    ext0 = 0

    # Taken from ACS. Should work for WFC3 too.
        
    x1_left  = 2027    ; x2_left  = 2048
    x1_right = x2_left ; x2_right = 2089
    xImageL = 2048; xImageR = xImageL * 2

    y1 = [0, 1548]; y2 = [500, 2048]  # near amp

    #-------------------------
    # Process image.

    pf = pyfits.open(outputfinal, mode='update')

    for i,ext in enumerate(ext_sci):

        #-------------------------
        # Extract strips near amp boundary.
    
        data_l = pf[ext].data[y1[i]:y2[i], x1_left:x2_left]
        data_r = pf[ext].data[y1[i]:y2[i], x1_right:x2_right]

        #-------------------------
        # Calculate and fit signal distributions.

        plotTitleL = '%s[%i] %s quad' % (outputfinal, ext, 'left')
        plotTitleR = '%s[%i] %s quad' % (outputfinal, ext, 'right')

        a = mysky.msky(data_l, do_plot=doplot, verbose=doplot, ptitle=plotTitleL)
        signalL = a[0]
        a = mysky.msky(data_r, do_plot=doplot, verbose=doplot, ptitle=plotTitleR)
        signalR = a[0]

        sDiff = signalL - signalR
        sDiff_abs = abs(sDiff)
        hJump = 0.5 * sDiff

        print '    %10s %10s %10s' % ('Left side', 'Right side', 'Half jump')
        print '    %10.6f %10.6f %10.6f' % (signalL, signalR, hJump)

        #-------------------------
        # Check threshold

        if sDiff_abs >= threshold:
            print '... removing the jump'
            pf[ext].data[:, 0:xImageL] -= hJump
            pf[ext].data[:, xImageL:xImageR] += hJump

            # Update header
            pf[ext0].header.add_history('Half jump %.6f corrected in EXT %i' % (hJump, ext))
        else:
            print '... no jump'

    pf.close()

def normalize_flat(outputfinal):
    """
    Rescales the level of the flat field to 1.

    This is not a requirement but has traditionally been
    done for the WFPC2 and ACS CCDs. A small box (101x101,
    like ACS) is defined in science extension 4 (the most
    used, we think) and its median is calculated. The four
    quadrants are then normalized to this median to
    preserve the relative offset in their levels (adjusted
    with the gains in calacs). The small box was defined
    in region free of droplets (by Elena Sabbi on Sep 23,
    2008). The exposure time keyword(s) are untouched.

    Extensions:
        * SCI (1 and 4) are normalized (each of the 4 quadrants)
        * ERR (2 and 5) are normalized (each of the 4 quadrants)
        * DQ (3 and 6) is unchanged

    Extension relationship with amplifiers:
        * [1], [2], [3] : C and D
        * [4], [5], [6] : A and B

    .. note:: `superflat` is disabled for ACS in `wfc_reference.py`. WFC3 needs to test if result is still the same with IRAF tasks replaced.

    Parameters
    ----------
    outputfinal: string
        Flatfield to normalize.

    """
    fd = pyfits.open(outputfinal, mode='update')

    # Get the binnings and image size.
    bin1 = fd[1].header['BINAXIS1']
    bin2 = fd[1].header['BINAXIS2']
    naxis1 = fd[1].header['NAXIS1'] # Horizontal/serial i.e. 4096
    naxis2 = fd[1].header['NAXIS2'] # Vertical/parallel i.e. 2051

    # Calculate the median in the box of extension 4 (amp A).
    if bin1 == 1 and bin2 == 1:
        im4med = fd[4].data[327:429,1031:1133]
    elif bin1 == 2 and bin2 == 2:
        im4med = fd[4].data[163:215,515:567]
    elif bin1 == 3 and bin2 == 3:
        im4med = fd[4].data[108:142,343:377]
        
    median_in_box = numpy.median( im4med )
    print '    Calculated a median of: %8.2f' % median_in_box

    # Normalize SCI and ERR by that median.
    for ext in (1,2,4,5): fd[ext].data /= median_in_box

    fd[0].header.add_history('Extensions 1,4 and 2,5 normalized with a level of: %8.2f DNs' % median_in_box)
    fd.close()

#------- Generic functions. -------

def before_combine(rawList, biasFile, refType, instrument, pctetab='', darkfile=''):
    """
    Adjusts keyword values in preparation for running
    `calwf3` or `calacs` to perform cosmic-ray rejection.

    These are applied to each frame in `rawList` to
    set calibration switches:
        * BIAS:
            * PERFORM: DQICORR, BLEVCORR, CRCORR
            * OMIT: all others
        * DARK:
            * PERFORM: DQICORR, BLEVCORR, BIASCORR, CRCORR
            * OMIT: all others
        * FLAT:
            * PERFORM: DQICORR, BLEVCORR, BIASCORR, DARKCORR, CRCORR
            * OMIT: all others

    If DQICORR is set to PERFORM, then the bad pixels
    flagged in BPIXTAB (for example, bad columns have a
    flag of 4) will have values of zero in the science
    extensions (1 and 4) of the output CRJ file (and
    consequently in the superbias).

    On the other hand, if DQICORR is set to OMIT, then
    the DQ extensions should be all set to zero and no
    pixel will be set to zero in the science extensions
    (1 and 4), so the images will look 'normal'.

    .. note:: `pctetab` is for ACS pixel-based CTE correction.

    Parameters
    ----------
    rawList: list
        List of RAW flat frames.
    
    biasFile: string
        Custom BIASFILE for BIASCORR.
    
    refType: {'bias', 'dark', 'flat'}
        Type of reference file in `rawList`.
    
    instrument: {'acs', 'wfc3'}
        Instrument name.

    pctetab: string
        Custom PCTETAB for PCTECORR. No CTE correction if
        not given. Use for DARK and FLAT, but not BIAS.

    darkfile: string
        Custom superdark for DARKCORR. Use for FLAT.
        If `pctetab` is also given, DRKCFILE is populated
        instead of DARKFILE.

    """    
    nFiles = len(rawList)
    hdrExt = 0

    # List of keywords to PERFORM
    perfList = ['DQICORR', 'BLEVCORR', 'CRCORR']

    # List of keywords to OMIT
    omitList = ['ATODCORR', 'FLSHCORR', 'EXPSCORR', 'SHADCORR', 'FLATCORR', 'PHOTCORR', 'DRIZCORR']

    # refType dependent keywords
    if refType == 'bias':
        omitList.append('BIASCORR')
        omitList.append('DARKCORR')
    elif refType == 'dark':
        perfList.append('BIASCORR')
        omitList.append('DARKCORR')
    else: # flat
        perfList.append('BIASCORR')
        perfList.append('DARKCORR')

    # instrument dependent keywords
    if instrument == 'acs':
        omitList.append('RPTCORR')
    
    for frame in rawList:
        fd = pyfits.open(frame, mode='update')
        fd[hdrExt].header.update('CRSPLIT', nFiles)

        # Custom superbias for BIASCORR
        if biasFile:
            fd[hdrExt].header.update('BIASFILE', biasFile)

        # Pixel-based CTE correction
        if pctetab:
            fd[hdrExt].header.update('PCTECORR', 'PERFORM')
            fd[hdrExt].header.update('PCTETAB', pctetab)
        else:
            fd[hdrExt].header.update('PCTECORR', 'OMIT')
 
        # Custom superdark for DARKCORR
        if darkfile:
            if pctetab:
                fd[hdrExt].header.update('DRKCFILE', darkfile)
            else:
                fd[hdrExt].header.update('DARKFILE', darkfile)

        # Set to PERFORM
        for hdrKey in perfList:
            fd[hdrExt].header.update(hdrKey, 'PERFORM')

        # Set to OMIT
        for hdrKey in omitList:
            fd[hdrExt].header.update(hdrKey, 'OMIT')

        fd.close()

def check_frames(rawList, instrument):
    """
    Print basic header keywords and set subarray flag
    to be returned.

    This subarray flag is set by checking NEXTEND of
    the last file in `rawList`.

    .. note:: WFC3 has to check whether its logic still holds.

    Parameters
    ----------
    rawList: list
        List file names for RAW frames.
    
    instrument: {'acs', 'wfc3'}
        Instrument name.

    Returns
    -------
    isSubArr: bool
        `True` if input frames are subarrays, else `False`.

    """
    print os.linesep, 'Will process these frames:', os.linesep
    hdr0, hdr1 = 0, 1
    comma = ','

    if instrument == 'wfc3':
        filKey = ['FILTER']
        fmt = '%-25s %6.1f %-8s %-4s %3.1f %ix%i %i %4i %4i'
    else: # acs
        filKey = ['FILTER1', 'FILTER2']
        fmt = '%-25s %6.1f %-17s %-4s %3.1f %ix%i %i %4i %4i'

    for frame in rawList:
        pf = pyfits.open(frame)
        
        filter = []
        for f in filKey: filter.append( pf[hdr0].header[f] )
        
        exptime  = pf[hdr0].header['EXPTIME']
        ccdamp   = pf[hdr0].header['CCDAMP']
        ccdgain  = pf[hdr0].header['CCDGAIN']
        nextend  = pf[hdr0].header['NEXTEND']
        bin1     = pf[hdr1].header['BINAXIS1']
        bin2     = pf[hdr1].header['BINAXIS2']
        naxis1   = pf[hdr1].header['NAXIS1']
        naxis2   = pf[hdr1].header['NAXIS2']

	print fmt % (frame, exptime, comma.join(filter), ccdamp, ccdgain, bin1, bin2, nextend, naxis1, naxis2)

        pf.close()

    # Set subarray flag. Assumes all frames in rawList have same format.
    # This logic only works for ACS/WFC or WFC3/UVIS.
    if nextend == 3:
        isSubArr = True
    else:       # 6
        isSubArr = False

    # Extra checks for WFC3 subarray
    if instrument == 'wfc3' and isSubArr:
        if naxis1 != 4142 or naxis2 != 2050 or ccdamp != 'B' or exptime != 0:
            raise Exception('This frame is not a UVIS1-2K4-SUB subarray.')

    return isSubArr

def buildasntable(rootname):
    """
    Build the association table.

    This table is used for `calwf3` or `calacs` in `combine`.

    Parameters
    ----------
    rootname: string
        Root of ASN table.

    Returns
    -------
    asn_file: string
        Filename of the generated ASN table.

    """
    asntab = asnutil.ASNTable(parseinput.parseinput('*raw.fits')[0], output=rootname)
    asntab.create()
    asntab.write()
    asn_file = rootname+'_asn.fits'

    # Change the MEMTYPE column (2nd column) to EXP-CRJ for each
    # exposure and the combined frame (last row) to PROD-CRJ.
    # (By default, buildasn sets these to EXP-DTH and PROD-DTH). 

    asn_table = pyfits.open(asn_file, mode='update')
    asn_table_data = asn_table[1].data

    nrows = asn_table_data.shape[0]
    nrows2 = nrows - 1

    for i in range(nrows2):
        asn_table_data.field('MEMTYPE')[i] = 'EXP-CRJ'
        product = asn_table_data.field('MEMNAME')[i]
        asn_table_data.field('MEMNAME')[i] = product

    asn_table_data.field('MEMTYPE')[nrows2] = 'PROD-CRJ'

    # Rewrite the product name in the last row to work around a bug in
    # PyFITS that leaves trailing blanks on the name (the process of
    # reading and rewriting the name strips off the blanks).
    # H.Bushouse, 22-June-2009
    product = asn_table_data.field('MEMNAME')[nrows2]
    asn_table_data.field('MEMNAME')[nrows2] = product

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # WFC3 UVIS had this commented out too
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #for i in range(nrows): asn_table_data.field('MEMPRSNT')[i] = 'no'

    asn_table.close()

    return asn_file

def combine(asn_filename, cal_func, save_temp):
    """
    Perform cosmic-ray rejection on an association.

    Parameters
    ----------
    asn_filename: string
        ASN table name.
    
    cal_func: {'calacs', 'calwf3'}
        Calibration task to use.
    
    save_temp: bool
        Keep temporary files?

    Returns
    -------
    crj, crc: string
        CRJ and CRC (if exists) output names. Raise
        exception if CRJ not found.

    """
    cal_exe = cal_func + '.e'
    cal_log = cal_func + '.log'
    cmd = [cal_exe, '-v']
    if save_temp: cmd.append('-s')
    cmd.append(asn_filename)
    flog = open(cal_log,'w')
    retval = subprocess.call(cmd, stdout=flog)
    flog.close()

    # Find CRJ
    # Python will raise exception if not found
    crj = glob.glob('*_crj.fits')[0]

    # Find CRC
    # Aug 2011 - Only made by CALACS if PCTECORR=PERFORM
    try:
        crc = glob.glob('*_crc.fits')[0]
    except IndexError:
        crc = None

    return crj, crc

def crx2out(crx, outputfinal, save_temp, doSubArr=False):
    """
    Prepare `combine` output for further processing.

    One or both of the DQ extensions of `crx` may be empty.
    If so, populate with all-zeros of the correct dimensions
    and data type.

    For CALACS 8.0 onwards, `crx` is in ELECTRONS (not DN)
    even when FLATCORR=OMIT. When BUNIT=ELECTRONS, convert
    data back to DN.

    Parameters
    ----------
    crx: string
        CRJ or CRC output from `combine`.

    outputfinal: string
        Name of modified `crx` to be the reference file.

    save_temp: bool
        Keep original `crx`?

    doSubArr: bool
        Use logic for subarray instead of full-frame?
        
    """
    if doSubArr: # Subarray
        dq_exts     = (3, )
    else:        # Full-frame
        dq_exts     = (3, 6)

    # Copy the cr-rejected image to outputfinal
    # Open outputfinal for edit
    shutil.copyfile(crx, outputfinal)
    fd = pyfits.open(outputfinal, mode='update')

    # Fix DQ extensions

    for i_dq in dq_exts:
        try:
            dat = pyfits.getdata(crx, i_dq)
        except:
            i_err = i_dq - 1

            fd[i_dq].data = (fd[i_err].data * 0).copy().astype('short')
            fd[i_dq].header.update('EXTNAME', 'DQ')

    # Convert electrons to DN.

    bunit = fd['SCI',1].header['BUNIT']

    if 'ELECTRON' in bunit:
        gain = get_gains(fd)
        apply_gain(fd, gain, 'div')

    # Close file and remove unwanted file
    fd.close()
    if not save_temp: os.remove(crx)

def add_overscans(outputfinal, template, instrument, doSubArr=False):
    """
    The CRJ file has its overscan trimmed off. This function
    re-inserts the overscan. The overscan regions are set to
    zeros. Regions are pre-defined.

    **Full Frame**

    * 1 or [sci,1] or CCDCHIP = 2 : C and D : overscan at the top
    * 2 or [err,1] or CCDCHIP = 2 : C and D : overscan at the top
    * 3 or [dq,1]  or CCDCHIP = 2 : C and D : overscan at the top
    * 4 or [sci,2] or CCDCHIP = 1 : A and B : overscan at the bottom
    * 5 or [err,2] or CCDCHIP = 1 : A and B : overscan at the bottom
    * 6 or [dq,2]  or CCDCHIP = 1 : A and B : overscan at the bottom

    *WFC3*

    The imaging areas are (see ISR 03-14, Bushouse):
        * 1x1 binning:
            * Amp A : [4][26:2073,20:2070]
            * Amp B : [4][2134:4181,20:2070]
            * Amp C : [1][26:2073,1:2051]
            * Amp D : [1][2134:4181,1:2051]  
        * 2x2 binning:
            * Amp A : [4][13:1035,10:1034]
            * Amp B : [4][1064:2086,10:1034]
            * Amp C : [1][13:1035,1:1025]
            * Amp D : [1][1064:2086,1:1025]
        * 3x3 binning:
            * Amp A : [4][9:690,7:689]
            * Amp B : [4][711:1392,7:689]
            * Amp C : [1][9:690,1:683]
            * Amp D : [1][711:1392,1:683]

    For 2x2 and 3x3 binnings, the CR-rejected output
    image from calwf3 includes the mixed rows/columns
    (ICD-47, Section 10.4 and ISR 03-14 do not include
    them; see Bushouse e-mails of Aug 20, 2008).

    The dimensions of the CR-rejected image are:
        * 1x1 : 4096x2051
        * 2x2 : 2048x1026
        * 3x3 : 1364x684

    *ACS*

    The imaging areas are (1x1 binning):
        * Amp A : [4][25:2072,21:2068]
        * Amp B : [4][2073:4120,21:2068]
        * Amp C : [1][25:2072,1:2048]
        * Amp D : [1][2073:4120,1:2048]

    **Subarray**

    * 1 or [sci,1]
    * 2 or [err,1]
    * 3 or [dq,1]

    *WFC3*

    The imaging areas for 1x1 binning are as below.
    The raw frames have 23 columns of overscans at
    each end. See ISR 2003-14, Bushouse.

    Amp B : [24:4119,1:2050] 

    *ACS*

    See Instrument Handbook for a list of supported
    subarrays and their respective layouts and
    dimensions.

    Parameters
    ----------
    outputfinal: string
        CR-rejected image to be processed. It is therefore
        identical to the original but with the overscans
        included (the idea is to leave the CR-rejected
        image as is, so it does not need to be regenerated).
    
    template: string
        One of the RAW files.
    
    instrument: {'acs', 'wfc3'}
        Instrument name.

    doSubArr: bool
        Use logic for subarray instead of full-frame?

    """
    ext0, ext1 = 0, 1
    fd = pyfits.open(outputfinal, mode='update')
    ft = pyfits.open(template)

    # Header info

    ccdamp   = fd[ext0].header['CCDAMP']
    n_ext    = fd[ext0].header['NEXTEND']

    binning  = fd[ext1].header['BINAXIS1']
    binning2 = fd[ext1].header['BINAXIS2']
    crj_xsz  = fd[ext1].header['NAXIS1']
    crj_ysz  = fd[ext1].header['NAXIS2']

    ccdchip  = ft[ext1].header['CCDCHIP']
    out_xsz  = ft[ext1].header['NAXIS1']
    out_ysz  = ft[ext1].header['NAXIS2']

    tpl_keys =['LTV1', 'LTV2', 'LTM1_1', 'LTM2_2']

    # Define the imaging areas.

    out_x1, out_x2, out_y1, out_y2 = {}, {}, {}, {}

    if instrument == 'wfc3':
        if doSubArr:
            if binning == 1:  # Only AMP B in wfc3_2k4sub_bias.py
                out_x1['B'], out_x2['B'], out_y1['B'], out_y2['B'] =   23, 4119,  0, 2050
            else:
                fd.close()
                ft.close()
                raise Exception("The subarray imaging areas are only defined for a binning of 1.")
        else:
            if binning == 1:
                out_x1['A'], out_x2['A'], out_y1['A'], out_y2['A'] =   25, 2073, 19, 2070
                out_x1['B'], out_x2['B'], out_y1['B'], out_y2['B'] = 2133, 4181, 19, 2070
                out_x1['C'], out_x2['C'], out_y1['C'], out_y2['C'] =   25, 2073,  0, 2051
                out_x1['D'], out_x2['D'], out_y1['D'], out_y2['D'] = 2133, 4181,  0, 2051
            elif binning == 2:  # Include the mixed rows/columns.
                out_x1['A'], out_x2['A'], out_y1['A'], out_y2['A'] =   13, 1037,  9, 1035
                out_x1['B'], out_x2['B'], out_y1['B'], out_y2['B'] = 1065, 2089,  9, 1035
                out_x1['C'], out_x2['C'], out_y1['C'], out_y2['C'] =   13, 1037,  0, 1026
                out_x1['D'], out_x2['D'], out_y1['D'], out_y2['D'] = 1065, 2089,  0, 1026
            elif binning == 3:  # Include the mixed rows/columns.
                out_x1['A'], out_x2['A'], out_y1['A'], out_y2['A'] =    9,  691,  6,  690
                out_x1['B'], out_x2['B'], out_y1['B'], out_y2['B'] =  711, 1393,  6,  690
                out_x1['C'], out_x2['C'], out_y1['C'], out_y2['C'] =    9,  691,  0,  684
                out_x1['D'], out_x2['D'], out_y1['D'], out_y2['D'] =  711, 1393,  0,  684
            else:
                fd.close()
                ft.close()
                raise Exception("The fullframe imaging areas are only defined for binnings of 1, 2, and 3.")
    else:            # acs
        if doSubArr:         # Subarray with 22 pix physical overscan
            if binning == 1:
                out_x1['A'], out_x2['A'], out_y1['A'], out_y2['A'] = 22, out_xsz, 0, out_ysz
                out_x1['B'], out_x2['B'], out_y1['B'], out_y2['B'] =  0, crj_xsz, 0, out_ysz
                out_x1['C'], out_x2['C'], out_y1['C'], out_y2['C'] = 22, out_xsz, 0, out_ysz
                out_x1['D'], out_x2['D'], out_y1['D'], out_y2['D'] =  0, crj_xsz, 0, out_ysz
            else:
                fd.close()
                ft.close()
                raise Exception("The subarray imaging areas are only defined for a binning of 1.")
        else:
            if binning == 1: # Normal in-flight full-frame
                out_x1['A'], out_x2['A'], out_y1['A'], out_y2['A'] =   24, 2072, 20, 2068
                out_x1['B'], out_x2['B'], out_y1['B'], out_y2['B'] = 2072, 4120, 20, 2068
                out_x1['C'], out_x2['C'], out_y1['C'], out_y2['C'] =   24, 2072,  0, 2048
                out_x1['D'], out_x2['D'], out_y1['D'], out_y2['D'] = 2072, 4120,  0, 2048
            else:
                fd.close()
                ft.close()
                raise Exception("The fullframe imaging areas are only defined for binnings of 1.")

    # Define the image areas for CRJ

    ext_list = range(1, n_ext+1)
    crj_x1, crj_x2, crj_y1, crj_y2 = {}, {}, {}, {}
    amp_sci = {}

    if doSubArr:
        amp_sci[ccdchip] = [ ccdamp ]
        crj_x1[ccdamp], crj_x2[ccdamp], crj_y1[ccdamp], crj_y2[ccdamp] = 0, crj_xsz, 0, crj_ysz
    else:
        amp_sci[1] = ['A', 'B']
        amp_sci[2] = ['C', 'D']
        xx2 = crj_xsz / 2
        crj_x1['A'], crj_x2['A'], crj_y1['A'], crj_y2['A'] =   0,     xx2, 0, crj_ysz
        crj_x1['B'], crj_x2['B'], crj_y1['B'], crj_y2['B'] = xx2, crj_xsz, 0, crj_ysz
        crj_x1['C'], crj_x2['C'], crj_y1['C'], crj_y2['C'] =   0,     xx2, 0, crj_ysz
        crj_x1['D'], crj_x2['D'], crj_y1['D'], crj_y2['D'] = xx2, crj_xsz, 0, crj_ysz

    # Create dummy array of new size and fill it with zeroes

    dummy_arr = numpy.zeros( (out_ysz, out_xsz) )

    # Add overscan to each quad of each ext

    for ext in ext_list:
        
        # Determine data type
        if fd[ext].header['BITPIX'] == 16:
            out_dtype = 'int16'
        else:
            out_dtype = 'float32'

        # Store data from _crj.fits
        crj_data = fd[ext].data.copy().astype(out_dtype)

        # Expand the dimension to include overscans

        fd[ext].data = dummy_arr.copy().astype(out_dtype)

        # Copy data back to pre-defined image areas

        try:
            cc = ft[ext].header['CCDCHIP']
        except KeyError:
            cc = cc_last
        cc_last = cc

        for amp in amp_sci[cc]: fd[ext].data[ out_y1[amp]:out_y2[amp], out_x1[amp]:out_x2[amp] ] = crj_data[ crj_y1[amp]:crj_y2[amp], crj_x1[amp]:crj_x2[amp] ].copy().astype(out_dtype)

        # Update some WCS header keywords from the template RAW image. calacs needs at least LTV1, LTV2, LTM1_1, and LTM2_2 to manage the appropriate overscan regions.
        for s in tpl_keys: fd[ext].header.update(s, ft[ext].header[s])

    fd.close()
    ft.close()

def get_gains(fd):
    """
    Get gains from CCDTAB.

    .. note:: Could add parse option for custom CCDTAB but not implemented.

    Parameters
    ----------
    fd: pointer
        FITS file pointer.

    Returns
    -------
    gain: dictionary
        Gain for each amp.

    """
    ext0, ext1 = 0, 1

    # Get header values for CCDTAB selection

    bin1 = fd[ext1].header['BINAXIS1']
    bin2 = fd[ext1].header['BINAXIS2']
    commanded_amps     = fd[ext0].header['CCDAMP']
    commanded_gain     = fd[ext0].header['CCDGAIN']
    commanded_ccdofsta = fd[ext0].header['CCDOFSTA'] # Assume other amps have the same offset.

    ccdtab = fd[ext0].header['CCDTAB']

    # Expand iref (WFC3) or jref (ACS)
    iref_dir = os.getenv('iref')
    jref_dir = os.getenv('jref')
    if 'iref$' in ccdtab:
        ccdtab = ccdtab.replace('iref$', iref_dir + os.sep)
    elif 'jref$' in ccdtab:
        ccdtab = ccdtab.replace('jref$', jref_dir + os.sep)

    if not os.path.isfile(ccdtab):
        fd.close()
        raise Exception('Missing CCDTAB file: %s.' % ccdtab)

    print '    Gains from:', ccdtab
    ftb = pyfits.open(ccdtab)
    ftb_data = ftb[ext1].data
    ftb_mask = ( (ftb_data.field('CCDAMP') == commanded_amps) & \
                 (ftb_data.field('CCDGAIN') == commanded_gain) & \
                 (ftb_data.field('CCDOFSTA') == commanded_ccdofsta) & \
                 (ftb_data.field('BINAXIS1') == bin1) & \
                 (ftb_data.field('BINAXIS2') == bin2) )
    selected_row = ftb_data[ftb_mask].copy()
    ftb.close()

    # Extract gain values from first matching row
    # (the rest are duplicates for gain)

    n_match = len( selected_row.field('ATODGNA') )
    if n_match < 1:
        fd.close()
        raise Exception('No match in CCDTAB (CCDAMP, CCDGAIN, CCDOFSTA, BINAXIS*).')

    gain = {}
    gain['A'] = selected_row.field('ATODGNA')[0]
    gain['B'] = selected_row.field('ATODGNB')[0]
    gain['C'] = selected_row.field('ATODGNC')[0]
    gain['D'] = selected_row.field('ATODGND')[0]
    print '    Found these gains: %4.2f  %4.2f  %4.2f  %4.2f' % \
        (gain['A'], gain['B'], gain['C'], gain['D'])

    return gain

def apply_gain(fd, gain, mul_flag):
    """
    Multiply or divide each quadrant by appropriate gain.
    This is done to both SCI and ERR.

    Logic depends on SUBARRAY:
        #. True  - Subarray logic.
        #. False - Full-frame logic.

    Parameters
    ----------
    fd: pointer
        FITS file pointer.

    gain: dictionary
        Gains from `get_gains`.

    mul_flag: {'mul','div'}
        Indicator to multiply or divide image by gain.
        Multiply converts DN to electrons.
        Divide converts electrons to DN.

    """
    ext0, ext1 = 0, 1

    n_ext    = fd[ext0].header['NEXTEND'] + 1
    isSubArr = fd[ext0].header['SUBARRAY']

    naxis1   = fd[ext1].header['NAXIS1']  # Horizontal/serial

    if mul_flag == 'mul':
        g = gain

        # Retained for historical reason
        for ext in range(n_ext): fd[ext].header.update('CCDGAIN', 1.0)

    elif mul_flag == 'div': # Invert gain
        g = {}
        for key,val in gain.iteritems(): g[key] = 1.0 / val

    else:
        fd.close()
        raise Exception('Apply gain failed. Can only multiply or divide.')

    # Assume no overscan. If it exists, will not be excluded.

    if isSubArr:
        print '    Nultiplying subarray by its gain.'

        ext_list = (1, 2) # Subarray SCI and ERR
        reg_amp  = ( (fd[ext0].header['CCDAMP'], ), \
                     (fd[ext0].header['CCDAMP'], ) )

        regx1 = (0, )
        regx2 = (naxis1, )

    else:
        print '    Multiplying each quadrant by its gain.'

        ext_list = (1, 2, 4, 5) # Full-frame SCI and ERR
        reg_amp = ( ('C', 'D'), ('C', 'D'), \
                    ('A', 'B'), ('A', 'B') )   # EXT 1 2 4 5

        half_frame = int(naxis1/2)
        regx1 = (0, half_frame)      # AorC, BorD
        regx2 = (half_frame, naxis1) # AorC, BorD


    regr = range(len(regx1))
    for i, ext in enumerate(ext_list):
        for j in regr:
            amp = reg_amp[i][j]
            fd[ext].data[:, regx1[j]:regx2[j]] *= g[amp]

def electrons_normalize(outputfinal, convert_e=1, norm_exptime=1):
    """
    Converts the multi-extension WFC/UVIS or ACS/WFC
    `outputfinal` from DN to e- and normalizes by the
    integration time.

    By default, both of these actions are taken. The input
    file can be any size (with or without overscans) but
    must have the multi-extension format.

    Operations are only done to SCI (1 and 4) and ERR
    (2 and 5). DQ (3 and 6) is unchanged.

    Extension relationship with amplifiers:
        * [1], [2], [3] : C and D
        * [4], [5], [6] : A and B

    Parameters
    ----------
    outputfinal: string
        Image to convert and normalize.

    conver_e: {0, 1}
        Convert DN to electrons?
            * 1: Yes
            * 0: No

    norm_exptime: {0, 1}
        Divide with EXPTIME?
            * 1: Yes
            * 0: No

    """
    ext0, ext1 = 0, 1
    ext_list = (1, 2, 4, 5) # Full-frame SCI and ERR

    if convert_e:
        print '  Converting to electrons.'
        fd = pyfits.open(outputfinal, mode='update')

        # Get gain for each quadrant
        gain = get_gains(fd)

        # Multiply each "half" of the extensions by the appropriate gain.
        apply_gain(fd, gain, 'mul')

	fd.close()

    if norm_exptime == 1:
        print '  Normalizing the SCI and ERR extensions (1, 2, 4, 5) by the integration time.'

        fd = pyfits.open(outputfinal, mode='update')
        exptime = fd[ext0].header['EXPTIME']

        for ext in ext_list: fd[ext].data /= exptime

        fd[ext0].header.update('EXPTIME', 1.0)
        fd[ext0].header.update('TEXPTIME', 1.0)
        fd.close()

def zero_DQ_arrays(outputfinal, doSubArr=False):
    """
    Zero the DQ extensions.

    For full-frame, these extensions are 3 and 6.
    For subarray, only 3.

    Parameters
    ----------
    outputfinal: string
        The final, cr-rejected reference file.

    doSubArr: bool
        Use logic for subarray instead of full-frame?

    """
    if doSubArr:
        ext_list = [3]
    else:
        ext_list = [3, 6]
    pf = pyfits.open(outputfinal, mode='update')
    for ext in ext_list: pf[ext].data *= 0
    pf.close()

def clean_header(outputfinal, instrument):
    """
    Delete, null, or change some header keywords in the
    final reference file.

    Parameters
    ----------
    outputfinal: string
        Image with header to clean.
    
    instrument: {'acs', 'wfc3'}
        Instrument name.
        
    """
    ext0 = 0
    fd = pyfits.open(outputfinal, mode='update')

    # Delete some keywords.
    # Will delete even if the keyword doesn't exist.

    keywords_to_delete = ['OSLAMBDA', 'OSBANDW', 'OSMRMODE', 'OSLAMP',
                          'OSXELMP', 'OSQTHLMP', 'OSFILT0', 'OSFILT1',
                          'OSFILT2', 'OSFSHTR', 'OSDSHTR', 'OSSWMIRR',
                          'OSPT', 'OSFIBRE', 'OSIMGPOS', 'OSPT_V1',
                          'OSPT_V2', 'OSPT_V3', 'OSPT_V2A',
                          'OSPT_V3A', 'OSOA_V2', 'OSOA_V3',
                          'OSOA_V2A', 'OSOA_V3A', 'OSPT_X', 'OSPT_Y',
                          'OSPT_Z', 'OSM3A1', 'OSM3A2', 'OSM4A1',
                          'OSM4A2', 'OSHENE', 'OSHENEDB', 'OSLD635',
                          'OSLD810', 'OSLD1064', 'OSLD1310',
                          'OSSTABLM', 'OSSTABLA', 'OSET', 'OSETFOC',
                          'OSFLUX', 'OSNDCORR', 'OSDETCTR',
                          'OSWVFRNT', 'OSCALFLE', 'OSSANITY',
                          'OSCOMMNT', 'OSSCRIPT', 'TVNUM', 'FLASHDUR',
                          'FLASHCUR', 'FLASHSTA', 'PATTERN1',
                          'P1_SHAPE', 'P1_PURPS', 'P1_NPTS',
                          'P1_PSPAC', 'P1_LSPAC', 'P1_ANGLE',
                          'P1_FRAME', 'P1_ORINT', 'PATTSTEP',
                          'DATE-OBS', 'COMMENT',
                          'SUNANGLE', 'MOONANGL', 'SUN_ALT',
                          'FGSLOCK', 'GYROMODE', 'REFFRAME', 'MTFLAG',
                          'PROPOSID', 'LINENUM', 'PR_INV_L', 'PR_INV_F', 'PR_INV_M',
                          'TIME-OBS', 'EXPEND',
                          'EXPFLAG', 'QUALCOM1', 'QUALCOM2', 'QUALCOM3',
                          'QUALITY', 'PA_V3', 'POSTARG1', 'POSTARG2',
                          'STATFLAG', 'WRTERR', 'DQICORR', 'ATODCORR',
                          'BLEVCORR', 'BIASCORR', 'FLSHCORR', 'CRCORR',
                          'EXPSCORR', 'SHADCORR', 'DARKCORR', 'FLATCORR',
                          'PHOTCORR', 'RPTCORR', 'DRIZCORR', 'PCTECORR',
                          'BPIXTAB', 'CCDTAB', 'ATODTAB', 'OSCNTAB',
                          'BIASFILE', 'FLSHFILE', 'CRREJTAB', 'SHADFILE',
                          'DARKFILE', 'PFLTFILE', 'DFLTFILE', 'LFLTFILE',
                          'PHOTTAB', 'GRAPHTAB', 'COMPTAB', 'IDCTAB',
                          'DGEOFILE', 'MDRIZTAB', 'CFLTFILE', 'SPOTTAB',
                          'T_SGSTAR', 'P1_CENTR', 'ASN_ID', 'ASN_TAB',
                          'RA_TARG', 'DEC_TARG',
                          'DRKCFILE', 'IMPHTTAB', 'D2IMFILE','NPOLFILE',
                          'SCALENSE', 'INITGUES', 'SKYSUB',
                          'SKYSUM', 'CRSIGMAS', 'CRRADIUS', 'CRTHRESH',
                          'BADINPDQ', 'REJ_RATE', 'CRMASK', 'MDRIZSKY']
    for keyword in keywords_to_delete: del fd[ext0].header[keyword]

    # Delete section titles in header

    comments_to_delete = ['      / POINTING INFORMATION',
                          '      / TARGET OFFSETS (POSTARGS)',
                          '      / OTFR KEYWORDS',
                          '      / PATTERN KEYWORDS',
                          '      / CALIBRATION REFERENCE FILES',
                          '      / CALIBRATION SWITCHES: PERFORM, OMIT, COMPLETE',
                          '      / ASSOCIATION KEYWORDS',
                          '      / PROPOSAL INFORMATION']    
    for i in reversed(range(len(fd[ext0].header))):
        if fd[ext0].header[i] in comments_to_delete:
            del fd[ext0].header[i+2]
            del fd[ext0].header[i+1]
            del fd[ext0].header[i]

    # Null some keywords.
    # Will add if the keyword doesn't exist.
    # Commented out as unnessecary.

    
    #keywords_to_null = ['PROPOSID', 'LINENUM', 'PR_INV_L', 'PR_INV_F', 'PR_INV_M',
    #                    'TIME-OBS', 'EXPSTART', 'EXPEND']
    keywords_to_null = ['EXPSTART'] 
    for keyword in keywords_to_null: fd[ext0].header.update(keyword, '')
    

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # THIS WAS ALREADY COMMENTED IN WFC3.
    # JUST LEAVE DQICORR AS COMPLETE, NOT OMIT.
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # Change some keywords.
    #fd[ext0].header.update('DQICORR', 'OMIT')  # If DQICORR = OMIT in the RAW files.

    fd.close()

def pad_dscrp(in_str):
    """
    Pad DESCRIP with dashes until required length is met.
  
    .. note:: The required length is 67 char, as stated by CDBS.

    Parameters
    ----------
    in_str: string
        String to pad.

    Returns
    -------
    out_str: string
        Padded string.

    """   
    dscrp_lgth = 67
    dash_char = '-'
    sz_str = len(in_str)
    
    if sz_str > dscrp_lgth:   # truncate
        out_str = in_str[0:dscrp_lgth]
    elif sz_str < dscrp_lgth: # pad
        out_str = in_str + dash_char*( dscrp_lgth - sz_str )
    else:                     # no change
        out_str = in_str

    return out_str

def get_annealstr(annealdate, instrument):
    """
    Construct DESCRIP string using CCD anneal info.
    
    .. note:: Only used by ACS but will work with WFC3. ACS is labeled as ACS-R after 2008.

    Parameters
    ----------
    annealdate: string
        Anneal date in the format of 'Mon DD YYYY HH:MM:SS'.

    instrument: {'acs', 'wfc3'}
        Instrument name.

    Returns
    -------
    out_str: string
        Constructed string.

    """
    out_str = ''

    # Define instrument label
    ins_str = instrument.upper()
    if instrument == 'acs':
        dval_in = time.strptime(annealdate, '%b %d %Y %H:%M:%S')
        if dval_in.tm_year > 2008: ins_str += '-R'

    # CCD anneal Mon DD YYYY HH:MM:SSUT XXX
    out_str = 'CCD anneal %sUT %s' % (annealdate, ins_str)

    return out_str

def cdbs_update_header(outputfinal, rawList, useafterdate, annealdate, refType, instrument, calName, userName, doSubArr=False):
    """
    Update the header of the reference file with appropriate
    CDBS keywords specified in TIR CDBS 2008-01.

    Parameters
    ----------
    outputfinal: string
        Image with header to update.
    
    rawList: list
        List of RAW files used to make `outputfinal`.

    useafterdate: string
        Optional value for USEAFTER. If None, will auto select
        from `rawList`. Format is 'YYYY-MM-DD' or
        'YYYY-MM-DD HH:MM:SS'.

    annealdate: string
        Anneal info for DESCRIP instead of normal description.
        Format is 'Mon DD YYYY HH:MM:SS'. DESCRIP will be 'CCD
        anneal Mon DD YYYY HH:MM:SSUT XXX'. For ACS, XXX is
        ACS-R if anneal date > 2008, else ACS. If None, will
        do normal description.

    refType: {'bias', 'dark', 'flat'}
        Type of reference file.
    
    instrument: {'acs', 'wfc3'}
        Instrument name.
    
    calName: {'calacs', 'calwf3'}
        Calibration task name.
    
    userName: string
        Name of reference file creator.

    doSubArr: bool
        Use logic for subarray instead of full-frame?

    """
    ext0, ext1 = 0, 1

    # Image extensions
    if doSubArr:
        ext_dat = [1, 2]
        ext_dq  = [3]
        if instrument == 'wfc3':
            txt_frame = '2K4, M512, and custom subarrays.'
        else: # acs
            txt_frame = 'subarrays in APERTURE.'
    else:
        ext_dat = [1, 2, 4, 5]
        ext_dq  = [3, 6]
        txt_frame = 'full-frame four amp readouts.'

    # Different date formats
    dfmt_ymd  = '%Y-%m-%d'
    dfmt_ymdt = '%Y-%m-%d %H:%M:%S'
    dfmt_dmy  = '%d/%m/%Y'
    dfmt_uaf  = '%b %d %Y %H:%M:%S'

    # For header update
    # after=500 (any large number can be used) is used to force keywords at the bottom of the header
    min_after_num = 500

    #-----------------------------------     
    # For the inflight pedigree, need the earliest and latest DATE-OBS
    # from the _raw.fits files.

    min_EXPSTART, max_EXPSTART = 999999.999, 0.0
    list_PROPOSID = []

    for frame in rawList:
        hdr_f0 = pyfits.getheader(frame, ext0)
        list_PROPOSID.append( hdr_f0['PROPOSID'] )
        mjd_f0 = hdr_f0['EXPSTART']

        if mjd_f0 < min_EXPSTART:
            min_EXPSTART = mjd_f0
            min_dateobs  = hdr_f0['DATE-OBS']
            min_timeobs  = hdr_f0['TIME-OBS']

        if mjd_f0 > max_EXPSTART:
            max_EXPSTART = mjd_f0
            max_dateobs  = hdr_f0['DATE-OBS']
            max_timeobs  = hdr_f0['TIME-OBS']

    # output format is: dd/mm/yyyy
    mindate = time.strftime(dfmt_dmy, time.strptime(min_dateobs, dfmt_ymd))  
    maxdate = time.strftime(dfmt_dmy, time.strptime(max_dateobs, dfmt_ymd))

    pedigree_inflight = '%s %s %s' % ('INFLIGHT', mindate, maxdate)

    #for History keyword with propID    
    uniq_PROPOSID = ''
    for pp in set(list_PROPOSID): uniq_PROPOSID += ' %i' % pp
    
    #For COMMENT
    if instrument == 'wfc3':
        comment_val = 'Reference file created by %s.' % userName
    else:  # acs
        comment_val = 'Super%s created by ACS autoreference pipeline' % (refType)





        

    #-----------------------------------     
    # For the inflight USEAFTER, update with date given on the command
    # line (-u switch) or default to mindate of the pedigree if no
    # date is provided.
    #
    # useafterdate can have time element, but not mindate.

    if useafterdate:
        if len(useafterdate.split()) > 1:  # YYYY-MM-DD HH:MM:SS
            ua_out = time.strptime(useafterdate, dfmt_ymdt)
        else:                              # YYYY-MM-DD
            ua_out = time.strptime(useafterdate, dfmt_ymd)
    else:
        if instrument == 'wfc3':
            ua_out = time.strptime(min_dateobs, dfmt_ymd)
        else:  # acs
            ua_out = time.strptime(min_dateobs+' '+min_timeobs, dfmt_ymdt)

    u_after = time.strftime(dfmt_uaf, ua_out)

    #-----------------------------------     
    # Now update the header of the output files.

    fd = pyfits.open(outputfinal, mode='update')
    
    del fd[ext0].header['HISTORY']

    #fd[ext0].header.add_history('The information listed above is output from %s' % calName)
    #fd[ext0].header.add_history('and indicates reference files used.')
    #fd[ext0].header.add_history(' ')
    fd[ext0].header.add_history(pad_dscrp(''))

    #-----------------------------------
    # Get pipeline version numbers

    cal_ver  = fd[ext0].header['CAL_VER']
    opus_ver = fd[ext0].header['OPUS_VER']

    # ----------
    # BIAS
    
    if refType == 'bias':
        after_num = 0

        fd[ext0].header.update('FILETYPE', 'BIAS')

        #---------------
        if annealdate:
            ds_str = get_annealstr(annealdate, instrument)
        else: # This is also true for ACS-R due to bias gradient in DSI.
            ds_str = 'Bias intended for use with %s' % txt_frame
        fd[ext0].header.update('DESCRIP', pad_dscrp(ds_str))

        #---------------     
        binning = fd[ext1].header['BINAXIS1']

        if instrument == 'wfc3':
            if doSubArr:
                fd[ext0].header.update('APERTURE', 'CHIP1_SUB_NOCORNERS')
            else:
                fd[ext0].header.update('CCDAMP', 'ABCD')
                if binning == 1:
                    fd[ext0].header.update('APERTURE', 'FULLFRAME_4AMP')
                else:
                    fd[ext0].header.update('APERTURE', 'ANY')
        else: # acs
            if not doSubArr:
                fd[ext0].header.update('CCDAMP', 'ABCDALL')  # not ANY
                fd[ext0].header.update('APERTURE', 'WFC')

        #---------------     
        for ext in ext_dat: # 1 2 4 5
            fd[ext].header.update('BUNIT', 'COUNTS')

    # ----------
    # DARK

    elif refType == 'dark':
        after_num = 0

        fd[ext0].header.update('FILETYPE', 'DARK') 
        fd[ext0].header.update('CCDAMP', 'ANY')

        #---------------
        if instrument == 'wfc3':
            ds_str = 'Dark created from in-flight WFC3/UVIS frames.'
            fd[ext0].header.update('SAMP_SEQ', 'NONE')
            fd[ext0].header.update('SUBTYPE', 'FullImag')
            fd[ext0].header.update('APERTURE', 'UVIS')
            for ext in ext_dat: # 1 2 4 5
                del fd[ext].header['BUNIT'] # Delete it because there's no e-/sec option.
        else: # acs
            ds_str = 'Dark created from in-flight ACS/WFC frames.'
            fd[ext0].header.update('APERTURE', 'WFC')
            fd[ext0].header.update('CCDGAIN', -1)
            for ext in ext_dat: # 1 2 4 5
                fd[ext].header.update('BUNIT', 'ELECTRONS/S')

        #---------------
        if annealdate: ds_str = get_annealstr(annealdate, instrument)
        fd[ext0].header.update('DESCRIP', pad_dscrp(ds_str))

    # ----------
    # FLAT (WFC3 only, disabled for ACS)
    
    elif refType == 'flat' and instrument == 'wfc3':
        after_num = min_after_num

        fd[ext0].header.update('APERTURE', 'UVIS')
        
        #---------------
        if annealdate:
            ds_str = get_annealstr(annealdate, instrument)
        else:
            ds_str = 'Flat created from in-flight WFC3/UVIS frames.'
        fd[ext0].header.update('DESCRIP', pad_dscrp(ds_str), after=after_num)

        #---------------
        if '_pfl.fits' in outputfinal:
            fd[ext0].header.update('FILETYPE', 'PIXEL-TO-PIXEL FLAT')
	elif '_dfl.fits' in outputfinal:
            fd[ext0].header.update('FILETYPE', 'DELTA FLAT')
	elif '_lfl.fits' in outputfinal:
            fd[ext0].header.update('FILETYPE', 'LARGE SCALE FLAT')
        else:
            fd.close()
            raise Exception('Unknown flatfield FILETYPE.')

        #---------------
        for ext in ext_dat: # 1 2 4 5
           fd[ext].header.update('BUNIT', 'UNITLESS')

    # ----------
    # UNKNOWN

    else:
        fd.close()
        raise Exception('No rules defined for reference type %s.' % refType)

    # Other common headers

    if after_num >= min_after_num:
        fd[ext0].header.update('PEDIGREE', pedigree_inflight, after=after_num)
        fd[ext0].header.update('USEAFTER', u_after, after=after_num)
    else:
        fd[ext0].header.update('PEDIGREE', pedigree_inflight)
        fd[ext0].header.update('USEAFTER', u_after)

    for ext in ext_dq: # 3 6
        fd[ext].header.update('BUNIT', 'UNITLESS')

    fd[ext0].header.add_history('This file was generated with wfc_reference.py' )
    fd[ext0].header.add_history('Raws taken from proposal%s' % (uniq_PROPOSID))
    fd[ext0].header.add_history('using %s version %s and %s' % (calName, cal_ver, opus_ver))
    fd[ext0].header.add_history('with the following files:')
    for frame in rawList: fd[ext0].header.add_history('%s' % frame) 

    # This adds COMMENT as in CDBS TIR 2009-01

    c_str = '= \''+comment_val+'\''
    fd[ext0].header.add_comment(c_str, before='origin')

    fd.close()
