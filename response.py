import ddosa
import astropy.io.fits as pyfits
import ogip.spec
from numpy import *
import sys
import os
import ddosa
import eddosa
import pilton
import dataanalysis.core as da
from dataanalysis.caches import cache_core
import timesystem as ts


#bash /Integral/throng/savchenk/projects/integral/ibismm/process/java/run_gen.sh response_3d_diff.fits  /Integral2/data/reduced/ddcache//byrev/1597/FinalizeLUT2.xvbase5.2//b4d67d69/lut2_1d_final.fits.gz  0.1


class ResponseRev(ddosa.DataAnalysis):
    version="v3"

    input_lut2=eddosa.FinalizeLUT2P4

    copy_cached_input=False

    cached=True

    r1=1
    r2=0
    depth_limit=0

    nchan=256
    
    def get_version(self):
        v=self.get_signature()+"."+self.version
        if self.nchan==62:
            v+="_62"
        v+="_r1_%.5lg"%self.r1+"_r2_%.5lg"%self.r2+("_d%.5lg"%self.depth_limit if self.depth_limit>0 else "")+".."
        return v


    def main(self):
        r1=self.r1
        r2=self.r2
        depth_limit=self.depth_limit

        cmd="bash "+os.environ['IBISMM_PROCESS_ROOT']+"/run_gen.sh "+os.environ['IBISMM_PROCESS_DATA']+"/response_3d_diff.fits "+self.input_lut2.lut2_1d.get_path()+" 0.1 %.5lg %.5lg %.5lg %i"%(r1,r2,depth_limit,self.nchan)
        print(cmd)
        os.system(cmd)

        self.out_response_3d_reconstructed=da.DataFile("out_response_3d_reconstructed.fits")
        self.out_response_3d_nores=da.DataFile("out_response_3d_nores.fits")
        self.out_response_3d_depth=da.DataFile("out_response_3d_depth.fits")
        self.out_response_3d=da.DataFile("out_response_3d.fits")

class OGIPResponse(ddosa.DataAnalysis):
    input_response=ResponseRev

    version="v2"

    cached=True
    copy_cached_input=False

    def main(self):
        input_file=self.input_response.out_response_3d_reconstructed.get_path()

        area_factor=450.

        if self.input_response.nchan==62:
            emin,emax=(lambda x:(x['E_MIN'],x['E_MAX']))(pyfits.open("/unsaved_data/savchenk/rmf_62bands.fits")['ISGR-EBDS-MOD'].data) # we should put this in IC or smth
        elif self.input_response.nchan==256:
            emin,emax=(lambda x:(x['E_MIN'],x['E_MAX']))(pyfits.open(os.environ['CURRENT_IC']+"/ic/ibis/mod/isgr_ebds_mod_0001.fits")['ISGR-EBDS-MOD'].data)

        print(("nchan",emin.shape))

        de=emax-emin
        #de=ones_like(de)*0.4787

        r3d=pyfits.open(input_file)[0].data #.transpose()[::-1,::-1]

        m_emin=arange(2466)*0.4787+13
        m_emax=(arange(2466)+1)*0.4787+13
        
 #       m_emin=arange(2200)*0.5+13
 #       m_emax=(arange(2200)+1)*0.5+13
        m_de=m_emax-m_emin

        def write_response(rt1, rt2):
            rmap=r3d[:,:,rt1:rt2].sum(2)/outer(m_de,ones_like(de))/area_factor
            #rmap=r3d[:,:,rt1:rt2].sum(2)/outer(de,ones_like(de))

            print((rmap.shape))

            rmap[:1,:]=0
            rmap[-2:,:]=0
            rmap[:,:1]=0
            rmap[:,-2:]=0

            key="rt%.5lg_%.5lg"%(rt1,rt2)

            rmffn="rmf_%s.fits"%key
            ogip.spec.RMF.from_arrays(m_emin,m_emax,rmap, emin, emax).to_fits(rmffn)

            #rmap_normalized=rmap/outer(ones_like(de),rmap.sum(1))

            #ogip.spec.RMF(emin,emax,emin,emax,rmap_normalized).write("rmf_normalized.fits")

            ogip.spec.PHAI.from_arrays(1, rate=rmap[50,:]*100,stat_err=rmap[50,:]*5).to_fits("pha_%s.fits"%key)
            return rmffn

        self.rmf_16_116=da.DataFile(write_response(16,116))
        self.rmf_16_50=da.DataFile(write_response(16,50))
        self.rmf_16_30=da.DataFile(write_response(16,30))

