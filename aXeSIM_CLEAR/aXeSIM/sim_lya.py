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
import axesim_wfc3_new
import scipy.interpolate
import pyfits
import matplotlib.pyplot as plt
import os
from astropy.cosmology import FlatLambdaCDM
from glob import glob
import scipy.optimize
import matplotlib


def simulate_many_noAGN():

    #orbits = [2,12, 40]
    redshift = [6.0, 6.5, 7.6]
    orbits = [12]
    #redshift = [7.0]
    #mag = [24]
    #EW = [50]
    #mag_str = ['24']
    mag = [24, 25, 26, 27]                    # J-band
    EW = [20, 50, 100, 150]                   # Figure 2, proposal, Angstroms
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
             
                    edit_onespec(z, cont_mag, kpc, pix)                                # Edit the onespec file for running aXeSIM
                    
                    w, d = spectrum_disk('disk.dat', cont_mag, EW_disk)                # Make the disk spectrum
               
                    file_ext = str(z)+'z_'+mag_str[i]+'mag_'+str(e)+'EW_'+str(orb)+'orb'
                    file_ext1 = file_ext.replace('.', 'p')
                    
                    #plt.plot(w, d) 
                    #plt.xlim(1180, 1260)
                    #plt.savefig('../input_specs/'+file_ext1+'input_spect')
                    #plt.close()
                    
                    if os.path.exists('../FITS_'+str(orb)+'orb/'+file_ext) == True:    # If the file exists already, delete it
                        os.system('rm -rf ../FITS_'+str(orb)+'orb/'+file_ext)
                        print 'Removed old version of ../FITS_'+str(orb)+'orb/'+file_ext
                    os.mkdir('../FITS_'+str(orb)+'orb/'+file_ext)                      # Naming scheme, probably shit
                    axesim_wfc3_new.axe(file_ext, orb)                                 # Run aXeSIM
                    interp_many(file_ext, orb)                                         # Interpolate the resulting spectrum
                    
    return


def edit_onespec(z, cont_mag, kpc, pix):      # spectemp is the line number of the fake galaxy in input_spectra.lis

    # read in the one_spec file
    fname = './save/one_spec_G102.lis'
    f = open(fname, 'w')
    f.write('# 1  NUMBER\n# 2  X_IMAGE\n# 3  Y_IMAGE\n# 4  A_IMAGE\n# 5  B_IMAGE\n# 6  THETA_IMAGE\n# 7  MAG_F1250W\n# 8  Z\n# 9  SPECTEMP\n# 10 MODIMAGE\n')

    x = [0, 240, 480, 720]
    y = [20, 60, 100, 140, 180, 220, 260, 300, 340, 380, 420, 460, 500, 540, 580, 620, 660, 700, 740, 780, 820, 860, 900, 940, 980]

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

            i += 1
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
    while i < 1900:
        j = 0
        p = 0
        while j < 1990:
            
            d = dnew[j:j+80, i:i+480]
            e = enew[j:j+80, i:i+480]

            if os.path.exists('../FITS_'+str(orb)+'orb/'+file_ext+'/driz_'+file_ext+'_'+str(k)+'.fits') == True: # If the file exists already, delete it
                os.remove('../FITS_'+str(orb)+'orb/'+file_ext+'/driz_'+file_ext+'_'+str(k)+'.fits')

            pyfits.writeto('../FITS_'+str(orb)+'orb/'+file_ext+'/driz_'+file_ext+'_'+str(k)+'.fits', (), header=g[0].header)
            pyfits.append('../FITS_'+str(orb)+'orb/'+file_ext+'/driz_'+file_ext+'_'+str(k)+'.fits', d)            
            pyfits.append('../FITS_'+str(orb)+'orb/'+file_ext+'/driz_'+file_ext+'_'+str(k)+'.fits', e)


            p += 1
            j += 80
            k += 1
        i += 480
            
    #os.remove('OUTSIM/'+file_ext+'_slitless_2.STP.fits')
    #os.remove('OUTSIM/'+file_ext+'_slitless_2.SPC.fits')
    #os.remove('OUTSIM/'+file_ext+'_images.fits')
    #os.remove('OUTSIM/'+file_ext+'_direct.fits')
    #os.remove('OUTSIM/'+file_ext+'_spectra.fits')
    #os.remove('OUTSIM/'+file_ext+'_slitless.fits')
    os.remove('DATA/'+file_ext+'_spectra.fits')
    os.remove('DATA/'+file_ext+'_images.fits')

    return


def gaussian(x,a,x0,sigma):
    
    return a*np.exp(-(x-x0)**2/(2*sigma**2))


def calc_Reff(mass, z):      # used to be van der Wel+ 2014 mass-size relation, currently just a norm. dist.

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
    cont = 10**(-(cont_mag+48.6)/2.5)
    # N.B. the wavelength below must be matched with the magnitude given in one_spec_G102.lis
    # i.e., currently it is set as MAG_F125W, so I use the pivot wavelength for that filter of 12486 A
    cont = cont * (3e18)/(12486**2)# * 1e18

    sigma = 2

    wave = np.arange(0, 2500, 0.5)
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


def line_fit(orb):

    x = open('../integrated_'+str(orb)+'orb.dat', 'w')
    x.write('#           ID                     Lya            err\n')
    files = glob('../FITS_'+str(orb)+'orb/*.dat')

    for file in files:

        wave, flux, err = np.loadtxt(file, unpack=True)
        flux = flux[np.logical_and(wave > 1130, wave < 1400)]
        wave = wave[np.logical_and(wave > 1130, wave < 1400)]
        err = err[np.logical_and(wave > 1130, wave < 1400)]
    
        flux_new = flux
        r = 0
        while r < 5:

            sort = sorted(flux_new)
            size = np.size(flux_new)
            bad = np.logical_and(flux_new > sort[int(0.15*size)], (flux_new < sort[int(0.85*size)]))
            fit = np.poly1d(np.polyfit(np.array(wave)[bad], np.array(flux_new)[bad], 1))
            #print len(flux_new), len(wave)
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

        #plt.step(wave, flux, 'b', lw = 1.5)
        #plt.plot(wave, fit_result, 'r-', linewidth = 2.5)
        #plt.xlim(1130, 1350)
        #plt.xlabel('Wavelength (A)')
        #plt.ylabel('Counts/sec')
        ##plt.ylabel('Flux (erg/cm$^2$/sec/A) x 10$^{-16}$')       # Maybe 10^-17?
        #plt.minorticks_on()
        #plt.savefig('../1Dspectra_'+str(orb)+'orb/'+name+'_fit.pdf', format='pdf')
        #plt.close()

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
    sens = pyfits.getdata('WFC3.IR.G102.1st.sens.2.fits')
    wavelength = sens['WAVELENGTH']
    sensit = sens['SENSITIVITY']
    
    for fol in folders:

        files = glob(fol+'/*')

        for j, f in enumerate(files):

            a = f.split('z')[0]
            z = float(a.split('/')[2])
            if z != 6.5:
                continue

            dat = pyfits.open(f)
            flux2d = dat[1].data
            err2d = dat[2].data
            ypix = len(flux2d[0])
            
            # Calculated by fitting a line to by-eye determination of 5 redshifts
            #c = -47.425*z + 6814.759
            wld = 24.5/2
            wl = 6710 + np.arange(0, (flux2d.shape)[1], 1)*wld - 15*wld
            wl0 = wl/(1+z)

            err2d[np.where(err2d == 0)] = 10*np.max(err2d)

            # Taken from Jon's IDL code - not the best way to do it I don't think
            #idx = []
            #wtest = wavelength/(1+z)
            #w = np.linspace(wtest[0], wtest[-1], 4751)
            #for i in wl0:
            #    a = np.abs(w - i).argmin()
            #    if a != a+2:
            #        idx.append(a)
            #se = sensit[idx]
            
            print wavelength/(1+z)
            print wl0
            func = scipy.interpolate.interp1d(wavelength/(1+z), sensit)
            se = func(wl0)

            flux2d_new = np.zeros(flux2d.shape)
            err2d_new = np.zeros(err2d.shape)
            for i, fl in enumerate(flux2d):
                if i == 0:
                    print fl/se
                flux2d_new[i] = fl
                err2d_new[i] = err2d[i]
                
            flux2d = flux2d_new
            err2d = err2d_new

            flux1d = []
            err1d = []
            for xx in np.arange(0, (flux2d.shape)[1], 1):
                flux1d.append(np.sum(flux2d[20:ypix-21,xx]))
                err1d.append(np.sum(err2d[20:ypix-21,xx]))

            plt.plot(wl0, flux1d)
            plt.show()
            
            ff = open(fol+'_'+str(j)+'.dat', 'w')
            ff.write('# Wavelength         Flux       Error\n')
            for i, el in enumerate(flux1d):
                ff.write(str(wl0[i])+'    '+str(el)+'   '+str(err1d[i])+'\n')

            ff.close()

    return            
            

def complete(orb):

    lya, error = np.loadtxt('../integrated_'+str(orb)+'orb.dat', unpack=True, usecols=(1,2))
    ID = np.genfromtxt('../integrated_'+str(orb)+'orb.dat', unpack=True, usecols=(0), dtype='str')

    name = []
    for i, n in enumerate(ID):
        temp = n.split('_'+str(orb)+'orb')
        name.append(temp[0])
    EE, PP, mm = [], [], []
    f = open('../completeness_'+str(orb)+'orb.dat', 'w')
    un = np.unique(name)  # Find the unique IDs in name
    for n in un:
        ind = [i for i, j in enumerate(name) if j == n]
        b = n.split('_')
        z = float(b[0].split('z')[0])
        mag = int(b[1].split('m')[0])
        ew = int(b[2].split('E')[0])
        flux = lya[ind]
        err = error[ind]
        #input_flux = (ew * 10**(-(mag-23.9)/2.5) * 10**(-29) * (3e18)/(4830**2)) # wrong wavelength for G102, btw
        
        a = 0        
        for j, fl in enumerate(flux):
            #if (fl > input_flux - 3*err[j]) and (fl < input_flux + 3*err[j]) and fl > 3*err[j]:
            if fl > 3*err[j]:
                a += 1.
        perc_det = (a/len(flux))
        print a
        f.write(n)
        f.write('    '+str(z))
        f.write('    '+str(mag))
        f.write('    '+str(ew))
        f.write('   '+str(perc_det)+'\n')
        #f.write('   '+str(np.median(r))+)
        #EE.append(ew)
        #PP.append(perc_det)
        #mm.append(mag)
    f.close()

    return# EE, PP, mm


def comp_hexplots(orb):

    z, mag, ew, comp = np.loadtxt('../completeness_'+str(orb)+'orb.dat', unpack=True, usecols=(1,2,3,4))

    #ind = np.where(ew == 200)
    
    fig, ax = plt.subplots()
    im = ax.hexbin(z, mag, gridsize = 3, C = comp, cmap = plt.cm.YlGnBu)
    cb = fig.colorbar(im, ax=ax)
    cb.set_label('Completeness')
    ax.set_xlabel(r'Redshift')
    ax.set_ylabel(r'Magnitude')
    ax.set_xlim(5.6, 8)
    ax.set_ylim(22.5, 28.5)
    im.set_clim(vmin =0.2,vmax=1)
    plt.savefig('../hex_plots/hex_z_mag_'+str(orb)+'orb_acc')
    plt.close()

    fig, ax = plt.subplots()
    im = ax.hexbin(ew, z, gridsize = 3, C = comp, cmap = plt.cm.YlGnBu)
    cb = fig.colorbar(im, ax=ax)
    cb.set_label('Completeness')
    ax.set_ylabel(r'Redshift')
    ax.set_xlabel(r'Restframe Equivalent Width')
    ax.set_ylim(5.3, 8.3)
    ax.set_xlim(-10, 180)
    im.set_clim(vmin =0.2,vmax=1)
    fig.savefig('../hex_plots/hex_z_EW_'+str(orb)+'orb_acc')
    plt.close(fig)

    fig, ax = plt.subplots()
    im = ax.hexbin(ew, mag, gridsize = 3, C = comp, cmap = plt.cm.YlGnBu)
    ax.set_xlabel(r'Restframe Equivalent Width')
    ax.set_ylabel(r'Magnitude')
    ax.set_xlim(-10, 180)
    ax.set_ylim(22.5, 28.5)
    cb = fig.colorbar(im, ax=ax)
    cb.set_label('Completeness')
    im.set_clim(vmin =0.2,vmax=1)
    fig.savefig('../hex_plots/hex_EW_mag_'+str(orb)+'orb_acc')
    plt.close(fig)

    return

