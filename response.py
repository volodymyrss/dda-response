import ddosa
import astropy.io.fits as fits
import heaspa
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

    input_lut2=eddosa.FinalizeLUT2

    copy_cached_input=False

    cached=True

    r1=1
    r2=0
    depth_limit=0

    
    def get_version(self):
        v=self.get_signature()+"."+self.version+"_r1_%.5lg"%self.r1+"_r2_%.5lg"%self.r2+("_d%.5lg"%self.depth_limit if self.depth_limit>0 else "")+".."
        return v


    def main(self):
        r1=self.r1
        r2=self.r2
        depth_limit=self.depth_limit

        cmd="bash "+os.environ['IBISMM_PROCESS_ROOT']+"/run_gen.sh response_3d_diff.fits "+self.input_lut2.lut2_1d.get_path()+" 0.1 %.5lg %.5lg %.5lg"%(r1,r2,depth_limit)
        print cmd
        os.system(cmd)

        self.out_response_3d_reconstructed=da.DataFile("out_response_3d_reconstructed.fits")
        self.out_response_3d_nores=da.DataFile("out_response_3d_nores.fits")
        self.out_response_3d_depth=da.DataFile("out_response_3d_depth.fits")
        self.out_response_3d=da.DataFile("out_response_3d.fits")

class OGIPResponse(ddosa.DataAnalysis):
    input_response=ResponseRev

    version="v2"

    cached=True

    def main(self):
        import heaspa
        input_file=self.input_response.out_response_3d_reconstructed.get_path()

        area_factor=450.

        emin,emax=(lambda
        x:(x['E_MIN'],x['E_MAX']))(pyfits.open(os.environ['CURRENT_IC']+"/ic/ibis/mod/isgr_ebds_mod_0001.fits")['ISGR-EBDS-MOD'].data)
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

        self.rmf_16_116=da.DataFile(write_response(16,116))
        self.rmf_16_50=da.DataFile(write_response(16,50))
        self.rmf_16_30=da.DataFile(write_response(16,30))


class ResponseIC(ddosa.DataAnalysis):
    write_caches=[cache_core.TransientCache,ddosa.MemCacheIntegralFallback]

    input_rev=ddosa.Revolution
    input_rsp=OGIPResponse

    cached=True

    ic_version=1

    ds_name="ISGR-RMF.-RSP"

    def fill_ds(self,f):
        d=pyfits.open(self.input_rsp.rmf_16_116.get_path())

        f[1].data=d[2].data
        #for k in ['DETCHANS']:
        for k in ['DETCHANS','NAXIS','NAXIS1','NAXIS2']:
            f[1].header[k]=d[2].header[k]
        f[1].header['TLMAX4']=256
        f[1].header['RISE_MAX']=116
        f[1].header['RISE_MIN']=16


    def main(self):
        out_fn=self.ds_name.lower().replace(".","").replace("-","_")+"_%.4i.fits"%int(self.input_rev.input_revid.str()) 

        dc=pilton.heatool("dal_create")
        dc["obj_name"]=out_fn
        dc["template"]=self.ds_name+".tpl"
        ddosa.remove_withtemplate(dc["obj_name"].value+"("+dc["template"].value+")")
        dc.run()

        val_start_ijd,val_stop_ijd=self.get_validity_ijd()

        val_start_utc=ts.converttime("IJD",val_start_ijd,"UTC")
        val_stop_utc=ts.converttime("IJD",val_stop_ijd,"UTC")

        f=pyfits.open(out_fn)

        self.fill_ds(f)


        f[1].header['ORIGIN']="Saclay/APC/ISDC"
        f[1].header['VERSION']=self.ic_version
        f[1].header['FILENAME']=out_fn
        f[1].header['LOCATN']=out_fn
        f[1].header['RESPONSI']="Volodymyr Savchenko"
        f[1].header['STRT_VAL']=val_start_utc
        f[1].header['END_VAL']=val_stop_utc
        f[1].header['VSTART']=val_start_ijd
        f[1].header['VSTOP']=val_stop_ijd
        f.writeto(out_fn,clobber=True)

        self.ic_ds_file=da.DataFile(out_fn)

        #import integralicindex as iii
        #ict=iii.ICTree()
        #ict.add(out_fn)
        #ict.write()


    def get_validity_ijd(self):
        return map(float,ts.converttime("REVNUM",self.input_rev.input_revid.str(),"IJD").split()[1:])



