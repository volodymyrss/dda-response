import ddosa
import pyfits
import heaspa
from numpy import *
import sys
import os
import ddosa
import eddosa
import dataanalysis as da

#bash /Integral/throng/savchenk/projects/integral/ibismm/process/java/run_gen.sh response_3d_diff.fits  /Integral2/data/reduced/ddcache//byrev/1597/FinalizeLUT2.xvbase5.2//b4d67d69/lut2_1d_final.fits.gz  0.1

class ResponseRev(ddosa.DataAnalysis):
    input_lut2=eddosa.FinalizeLUT2

    copy_cached_input=False

    cached=True

    def main(self):
        os.system("bash /Integral/throng/savchenk/projects/integral/ibismm/process/java/run_gen.sh response_3d_diff.fits "+self.input_lut2.lut2_1d.get_path()+" 0.1")

        self.out_response_3d_reconstructed=da.DataFile("out_response_3d_reconstructed.fits")

class OGIPResponse(ddosa.DataAnalysis):
    input_response=ResponseRev

    cached=True

    def main(self):
        input_file=self.input_response.out_response_3d_reconstructed.get_path()

        area_factor=450.

        #emin,emax=(lambda x:(x['E_MIN'],x['E_MAX']))(pyfits.open("/Integral/data/resources/rmf_62bands.fits")[3].data)
        emin,emax=(lambda x:(x['E_MIN'],x['E_MAX']))(pyfits.open("/Integral/data/resources/rmf_256bins.fits")['EBOUNDS'].data)
        de=emax-emin
        #de=ones_like(de)*0.4787

        r3d=pyfits.open(input_file)[0].data #.transpose()[::-1,::-1]

       # m_emin=arange(2466)*0.4787+13
       # m_emax=(arange(2466)+1)*0.4787+13
        
        m_emin=arange(2200)*0.5+13
        m_emax=(arange(2200)+1)*0.5+13
        m_de=m_emax-m_emin

        def write_response(rt1, rt2):
            rmap=r3d[:,:,rt1:rt2].sum(2)/outer(m_de,ones_like(de))/area_factor
            #rmap=r3d[:,:,rt1:rt2].sum(2)/outer(de,ones_like(de))

            print rmap.shape

            rmap[:1,:]=0
            rmap[-2:,:]=0
            rmap[:,:1]=0
            rmap[:,-2:]=0

            key="rt%.5lg_%.5lg"%(rt1,rt2)

            rmffn="rmf_%s.fits"%key
            heaspa.RMF(emin,emax,m_emin,m_emax,rmap).write(rmffn)

            #rmap_normalized=rmap/outer(ones_like(de),rmap.sum(1))

            #heaspa.RMF(emin,emax,emin,emax,rmap_normalized).write("rmf_normalized.fits")

            heaspa.PHA(rmap[50,:]*100,rmap[50,:]*5,1).write("pha_%s.fits"%key)
            return rmffn

        write_response(16,93)
        self.rmf_16_116=da.DataFile(write_response(16,116))
        write_response(16,50)
        write_response(50,116)
        write_response(0,255)
