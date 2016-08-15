###########################################################################################
## This script and all it's functions should run all the simulations through aXeSIM      ##
##                                                                                       ##
## Joanna Bridge, 2015                                                                   ##
## Revamped 2016 for Ly-alpha line simulation                                            ##
##                                                                                       ##
## Notes: 1.) All runs of aXeSIM use the 'exponential_disk.fits' and the 1400nm PSF      ##
##            fits files, i.e., no need to use mkobjects code                            ##
##                                                                                       ##
##                                                                                       ##
###########################################################################################

import numpy as np
#import axesim_wfc3_new
import scipy.interpolate
import pyfits
import matplotlib.pyplot as plt
import os
from astropy.cosmology import FlatLambdaCDM
from glob import glob
import scipy.optimize
import matplotlib
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
import matplotlib.gridspec as gridspec


def simulate_many_noAGN():

    orbits = [2,12, 40]
    redshift = [5.6, 6.0, 6.5, 7.0, 7.6]
    mag = [24, 25, 26, 27]                    # J-band
    EW = [20, 50, 100, 150, 200]     # Figure 2, proposal, Angstroms
    mag_str = ['24', '25', '26', '27']
    
    for z in redshift:
        for i, m in enumerate(mag):
            for j, e in enumerate(EW):
                for l, orb in enumerate(orbits):
                
                    cont_mag = m                  
                    EW_disk = e             

                    kpc, pix = [],[]
                    for b in xrange(100):
                        d, g = calc_Reff(m, z)
                        kpc.append(d)
                        pix.append(g)
             
                    edit_onespec_noAGN(z, cont_mag, kpc, pix)                          # Edit the onespec file for running aXeSIM
                    
                    w, d = spectrum_disk('disk.dat', cont_mag, EW_disk)                # Make the disk spectrum
               
                    file_ext = str(z)+'z_'+mag_str[i]+'mag_'+str(e)+'EW_'+str(orb)+'orb'
                    file_ext1 = file_ext.replace('.', 'p')
                    
                    plt.plot(w, d)
                    plt.xlim(1180, 1260)
                    plt.savefig('../input_specs/'+file_ext1+'input_spect')
                    plt.close()
                    
                    #if os.path.exists('../FITS_'+str(orb)+'orb/'+file_ext) == True:    # If the file exists already, delete it
                    #    os.system('rm -rf ../FITS_'+str(orb)+'orb/'+file_ext)
                    #    print 'Removed old version of ../FITS_'+str(orb)+'orb/'+file_ext
                    #os.mkdir('../FITS_'+str(orb)+'orb/'+file_ext)                      # Naming scheme, probably shit
                    #axesim_wfc3_new.axe(file_ext, orb)                                 # Run aXeSIM
                    #interp_many(file_ext, orb)                                         # Interpolate the resulting spectrum
                    
    return


def edit_onespec_noAGN(z, cont_mag, kpc, pix):      # spectemp is the line number of the fake galaxy in input_spectra.lis

    # read in the one_spec file
    fname = './save/one_spec_G102.lis'
    f = open(fname, 'w')
    f.write('# 1  NUMBER\n# 2  X_IMAGE\n# 3  Y_IMAGE\n# 4  A_IMAGE\n# 5  B_IMAGE\n# 6  THETA_IMAGE\n# 7  MAG_F475W\n# 8  Z\n# 9  SPECTEMP\n# 10 MODIMAGE\n')

    x = [0, 200, 400, 600, 800]
    y = [20, 60, 100, 140, 180, 220, 260, 300, 340, 380, 420, 460, 500, 540, 580, 620, 660, 700, 740, 780]

    i=0
    n=0
    for k in x:
        for u, l in enumerate(y):
            
            blah = np.arange(0.5, 10, 0.1)
            ind = [j for j, a in enumerate(blah) if "{0:.1f}".format(a) == "{0:.1f}".format(kpc[n])]
            if not ind:
                ind = [0]
                pix[n] = 0.5
            
            f.write(str(i+1))
            f.write(' '+str(k)+' '+str(l))
            f.write(' '+"{0:.2f}".format(2*pix[n])+' '+"{0:.2f}".format(2*pix[n]))
            f.write(' 90.0 ')
            f.write(str(cont_mag))
            f.write(' '+str(z)+' ')
            f.write('1 '+str(2+ind[0])+'\n')

            i += 2
            n += 1
    return 


def interp_many(file_ext, orb):   #Drizzle fits file from aXeSIM
    
    g = pyfits.open('OUTSIM/'+file_ext+'_slitless.fits')
    dat = g[1].data
    err = g[2].data
    
    x = np.arange(0, dat.shape[1])
    y = np.arange(0, dat.shape[0])
    xnew = np.arange(0, dat.shape[1]*2)/2.
    ynew = np.arange(0, dat.shape[0]*2)/2.
    
    f = scipy.interpolate.interp2d(x, y, dat, kind = 'cubic')
    dnew = f(xnew, ynew)
    dnew = (dnew * 90220 - np.median(dnew))/90220
    
    e = scipy.interpolate.interp2d(x, y, err, kind = 'cubic')
    enew = e(xnew, ynew)
    
    i, k = 0, 1
    while i < 2000:
        j = 0
        p = 0
        while j < 1540:
            
            d = dnew[j:j+82, i:i+400]
            e = enew[j:j+82, i:i+400]

            #if os.path.exists('../FITS_'+str(orb)+'orb/'+file_ext+'/driz_'+file_ext+'_'+str(k)+'.fits') == True: # If the file exists already, delete it
            #    os.remove('../FITS_'+str(orb)+'orb/'+file_ext+'/driz_'+file_ext+'_'+str(k)+'.fits')

            #pyfits.writeto('../FITS_'+str(orb)+'orb/'+file_ext+'/driz_'+file_ext+'_'+str(k)+'.fits', (), header=g[0].header)
            #pyfits.append('../FITS_'+str(orb)+'orb/'+file_ext+'/driz_'+file_ext+'_'+str(k)+'.fits', d)            
            #pyfits.append('../FITS_'+str(orb)+'orb/'+file_ext+'/driz_'+file_ext+'_'+str(k)+'.fits', e)
            pyfits.writeto('../'+file_ext+'_'+str(i)+str(k)+'.fits', (), header=g[0].header)
            pyfits.append('../'+file_ext+'_'+str(i)+str(k)+'.fits', d)            
            #pyfits.append('../FITS_'+str(orb)+'orb/'+file_ext+'/driz_'+file_ext+'_'+str(k)+'.fits', e)

            p += 1
            j += 80
            k += 1
        i += 400
            
    os.remove('OUTSIM/'+file_ext+'_slitless_2.STP.fits')
    os.remove('OUTSIM/'+file_ext+'_slitless_2.SPC.fits')
    os.remove('OUTSIM/'+file_ext+'_images.fits')
    os.remove('OUTSIM/'+file_ext+'_direct.fits')
    #os.remove('OUTSIM/'+file_ext+'_spectra.fits')
    #os.remove('OUTSIM/'+file_ext+'_slitless.fits')
    os.remove('DATA/'+file_ext+'_spectra.fits')
    os.remove('DATA/'+file_ext+'_images.fits')

    return


def gaussian(x,a,x0,sigma):
    
    return a*np.exp(-(x-x0)**2/(2*sigma**2))


def calc_Reff(mass, z):      # van der Wel+ 2014 mass-size relation                                                

    kpc = np.random.normal(2.25, 1.5)
    if z == 5.6:
        b = 5.927
    if z ==  6.0:
        b = 5.710
    if z == 6.5:
        b = 5.458
    if z == 7:
        b = 5.226
    if z == 7.6:
        b = 4.972
        
    arcs = kpc/b
    pix = arcs/0.12825       # 0.12825 "/pix 

    return kpc, pix


def spectrum_disk(filename, cont_mag, EW):    # Makes gaussians

    # convert magnitude to f_lambda continuum
    cont = 10**(-(cont_mag-23.9)/2.5) * 10**(-29)
    # N.B. the wavelength below must be matched with the magnitude given in one_spec_G141.lis
    # i.e., currently it is set as MAG_F475W, so I use the pivot wavelength for that filter of 4830 A
    cont = cont * (3e18)/(4830**2)# * 1e18

    sigma = 2

    wave = np.arange(0, 2000, 0.5)
    a = 1/np.sqrt(2*np.pi*sigma**2)
    g = gaussian(wave, a, 1216, sigma)   
    C = EW * cont
    
    Lya = C * g
    G = Lya + cont
    
    f = open('SIMDATA/'+filename, 'w')
    for i, w in enumerate(wave):
        f.write('%.5e  ' % w)
        f.write('%.5e\n' % G[i])
    f.close()

    return wave, G


def edit_onespec(z, cont_mag):      # spectemp is the line number of the fake galaxy in input_spectra.lis
    
    # read in the one_spec file
    fname = './save/one_spec_G102.lis'
    f = open(fname, 'r')
    lines = f.readlines()
    
    # replace the data with what I choose
    line = lines[-1].split(' ')
    line[7] = str(z)
    lines[-1] = ' '.join(line)

    line = lines[-2].split(' ')
    line[6] = str(cont_mag)
    line[7] = str(z)
    lines[-2] = ' '.join(line)
    
    f.close()

    # rewrite the lines to the same file
    f = open(fname, 'w')
    for l in lines:
        f.write(l)

    f.close()

    return 


def line_fit(orb):

    x = open('../integrated_'+str(orb)+'orb.dat', 'w')
    x.write('#           ID                     Lya            err\n')
    files = glob('../FITS_'+str(orb)+'orb/*.dat')

    for file in files:

        wave, flux, err = np.loadtxt(file, unpack=True)
        flux = flux[np.logical_and(wave > 1000, wave < 1400)]
        wave = wave[np.logical_and(wave > 1000, wave < 1400)]
        err = err[np.logical_and(wave > 1000, wave < 1400)]
    
        flux_new = flux
        r = 0
        while r < 10:

            sort = sorted(flux_new)
            size = np.size(flux_new)
            bad = np.logical_and(flux_new > sort[int(0.15*size)], (flux_new < sort[int(0.85*size)]))
            fit = np.poly1d(np.polyfit(np.array(wave)[bad], np.array(flux_new)[bad], 1))
            flux_new = flux_new - fit(wave)
            r += 1
        flux = flux_new

        global N
        N = 1
        global line
        line = np.array([1216])

        name = file.strip('../FITS_'+str(orb)+'orb/')
        
        p0 = init(line, flux)

        try:
            coeff, var_mat = scipy.optimize.curve_fit(func, wave, flux, p0 = p0, sigma=err, absolute_sigma=True)
        except RuntimeError:
            print name
            x.write(file+'  \n')
            continue

        mean_Lya = coeff[N+1]
        sig = coeff[-1]
        err = np.sqrt(np.diag(var_mat))
        err_Lya = err[N+1]
    
        fit_Lya = gaussian(wave, mean_Lya, line[0], sig)
        area_Lya = (np.trapz(fit_Lya, wave))
        fit_result = func(wave, *coeff)

        x.write(name+'  ')
        x.write(str(area_Lya)+'   ')
        x.write(str(err_Lya)+'\n')

        plt.step(wave, flux, 'b', lw = 1.5)
        plt.plot(wave, fit_result, 'r-', linewidth = 2.5)
        plt.xlim(1100, 1350)
        plt.xlabel('Wavelength (A)')
        plt.ylabel('Flux (erg/cm$^2$/sec/A) x 10$^{-16}$')       # Maybe 10^-17?
        plt.minorticks_on()
        plt.savefig('../1Dspectra_'+str(orb)+'orb/'+name+'_fit.pdf', format='pdf')
        plt.close()

    x.close()

    return


def func(x, *p):

    n = N

    poly = 0 * x
    for i in xrange(n + 1):
        poly = poly + p[i] * x**i

    gauss1 = p[1+n]*np.exp(-(x-line[0])**2/(2*(p[-1]**2)))      # Lya
    
    return gauss1 + poly


def init(line, flux):
    
    n = N
    p0 = np.zeros(shape = (n + len(line) + 2))
    p0[0] = np.mean(flux)
    p0[-1] = 50
    
    return p0


def distribution(mean, err):
    
    return np.random.normal(mean, err, 1000)


def spect1D(orb):

    folders = glob('../FITS_'+str(orb)+'orb/*')
    sens = pyfits.getdata('WFC3.IR.G102.1st.sens.1.fits')
    wavelength = sens['WAVELENGTH']
    sensit = sens['SENSITIVITY']
    fun = scipy.interpolate.interp1d(wavelength, sensit)
    
    for fol in folders:

        files = glob(fol+'/*')

        for j, f in enumerate(files):

            a = f.split('z')[0]
            z = float(a.split('/')[2])

            dat = pyfits.open(f)
            flux2d = dat[1].data
            err2d = dat[2].data
            ypix = len(flux2d[0])
            wld = 46.5/2

            # Calculated by fitting a line to by-eye determination of 5 redshifts
            c = 1222*z + 684.58200637
            wl = c + np.arange(0, (flux2d.shape)[0], 1)*wld - 15*wld
            wl0 = wl/(1+z)

            err2d[np.where(err2d == 0)] = 10*np.max(err2d)
            print flux2d.shape

            idx = np.abs(wl - wavelength[0]).argmin()
            print idx
            wl = wl[idx:]
            flux2d = flux2d[:, idx:]
            err2d = err2d[:, idx:]
            
            se = fun(wl)

            flux1d = []
            err1d = []
            for xx in np.arange(0, (flux2d.shape)[0], 1):
                flux1d.append(np.sum(flux2d[xx,8:ypix-9]))
                err1d.append(np.sum(err2d[xx,8:ypix-9]))
            print np.array(flux1d).shape
            ff = open(fol+'_'+str(j)+'.dat', 'w')
            ff.write('# Wavelength         Flux       Error\n')
            for i, el in enumerate(flux1d):
                if i > len(se)-1:
                    continue
                else:
                    ff.write(str(wl0[i])+'    '+str(el/se[i])+'   '+str(err1d[i]/se[i])+'\n')
                    #ff.write(str(wl0[i])+'    '+str(el)+'   '+str(err1d[i])+'\n')

            ff.close()

    return            
            

def complete(orb):

    lya, error = np.loadtxt('../integrated_'+str(orb)+'orb.dat', unpack=True, usecols=(1,2))
    ID = np.genfromtxt('../integrated_'+str(orb)+'orb.dat', unpack=True, usecols=(0), dtype='str')

    name = []
    for i, n in enumerate(ID):
        temp = n.split('_'+str(orb)+'orb')
        name.append(temp[0])
    
    f = open('../acc_'+str(orb)+'orb.dat', 'w')
    un = np.unique(name)  # Find the unique IDs in name
    for n in un:
        ind = [i for i, j in enumerate(name) if j == n]
        b = n.split('_')
        z = float(b[0].split('z')[0])
        mag = int(b[1].split('m')[0])
        ew = int(b[2].split('E')[0])
        flux = lya[ind]
        err = error[ind]
        input_flux = (ew * 10**(-(mag-23.9)/2.5) * 10**(-29) * (3e18)/(4830**2))
        
        a = 0        
        for j, fl in enumerate(flux):
            if (fl > input_flux - 3*err[j]) and (fl < input_flux + 3*err[j]) and fl > 3*err[j]:
                a += 1.
                print fl, input_flux, err[j]
        perc_det = (a/len(flux))
        print a
        f.write(n)
        f.write('    '+str(z))
        f.write('    '+str(mag))
        f.write('    '+str(ew))
        f.write('   '+str(perc_det)+'\n')
        #f.write('   '+str(np.median(r))+)

    f.close()

    return 


def comp_hexplots(orb):

    z, mag, ew, comp = np.loadtxt('../completeness_'+str(orb)+'orb.dat', unpack=True, usecols=(1,2,3,4))

    ind = np.where(ew == 200)
    
    fig, ax = plt.subplots()
    im = ax.hexbin(z[ind], mag[ind], gridsize = 4, C = comp[ind], cmap = plt.cm.YlGnBu)
    cb = fig.colorbar(im, ax=ax)
    cb.set_label('Completeness')
    ax.set_xlabel(r'Redshift')
    ax.set_ylabel(r'Magnitude')
    ax.set_xlim(5.2, 8)
    ax.set_ylim(23, 28)
    im.set_clim(vmin =0,vmax=1)
    plt.savefig('../hex_plots/hex_z_mag_'+str(orb)+'orb_acc')
    plt.close()

    fig, ax = plt.subplots()
    im = ax.hexbin(ew, z, gridsize = 4, C = comp, cmap = plt.cm.YlGnBu)
    cb = fig.colorbar(im, ax=ax)
    cb.set_label('Completeness')
    ax.set_ylabel(r'Redshift')
    ax.set_xlabel(r'Equivalent Width')
    ax.set_ylim(5.2, 8)
    ax.set_xlim(25, 225)
    im.set_clim(vmin =0,vmax=1)
    fig.savefig('../hex_plots/hex_z_EW_'+str(orb)+'orb_acc')
    plt.close(fig)

    fig, ax = plt.subplots()
    im = ax.hexbin(ew, mag, gridsize = 4, C = comp, cmap = plt.cm.YlGnBu)
    ax.set_xlabel(r'Equivalent Width')
    ax.set_ylabel(r'Magnitude')
    ax.set_xlim(25, 225)
    ax.set_ylim(23, 28)
    cb = fig.colorbar(im, ax=ax)
    cb.set_label('Completeness')
    im.set_clim(vmin =0,vmax=1)
    fig.savefig('../hex_plots/hex_EW_mag_'+str(orb)+'orb_acc')
    plt.close(fig)

    return

