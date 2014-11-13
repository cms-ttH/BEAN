<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">
<html>
<META HTTP-EQUIV="CACHE-CONTROL" CONTENT="NO-CACHE">
<META HTTP-EQUIV="PRAGMA" CONTENT="NO-CACHE">
<META HTTP-EQUIV="EXPIRES" CONTENT="0">
<head>
<title>CMSSW/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py</title>
<base href="https://cmssdt.cern.ch/SDT/lxr/">
<link href="lxr.css" rel="STYLESHEET" type="text/css">
<script type="text/javascript">
function ensureFocus()
{
  if (document.getElementById("focus"))
	  {
	  document.getElementById("focus").focus();
	  }
}
</script>
</head>


<body onload="ensureFocus()">

<table width='100%' border='0' cellpadding='0' cellspacing='0'>
  <tr>
    <td valign='top'>
      <!-- put local logo or links here -->
    </td>
    <td>
      <table width='100%' border='0' cellpadding='0' cellspacing='0'>
        <tr>
          <td align='center'>
            <h1>The LXR Cross Referencer</h1>
          </td>
        </tr>
        <tr>
          <td align="center"><span class="banner"><a class='banner' href="/SDT/lxr/source/?v=CMSSW_5_3_21">CMSSW</a>/ <a class='banner' href="/SDT/lxr/source/SimGeneral/?v=CMSSW_5_3_21">SimGeneral</a>/ <a class='banner' href="/SDT/lxr/source/SimGeneral/MixingModule/?v=CMSSW_5_3_21">MixingModule</a>/ <a class='banner' href="/SDT/lxr/source/SimGeneral/MixingModule/python/?v=CMSSW_5_3_21">python</a>/ <a class='banner' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21">mix_2012A_Profile_PoissonOOTPU_cfi.py</a></span></td>
        </tr>
      </table>
    </td>
    <td align='right'>
      
      [&nbsp;<span class='modes-sel'>source navigation</span>&nbsp;]<br>
      [&nbsp;<a class='modes' href="/SDT/lxr/diff/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21">diff markup</a>&nbsp;]<br>
      [&nbsp;<a class='modes' href="/SDT/lxr/ident?v=CMSSW_5_3_21">identifier search</a>&nbsp;]<br>
      [&nbsp;<a class="modes" href="/SDT/lxr/search?v=CMSSW_5_3_21">general search</a>&nbsp;]<br>
    </td>
  </tr>
  <tr><td colspan='3'>&nbsp;</td></tr>
  <tr>
   <td>&nbsp;</td>
   <td colspan='2'>
     <table width="100%" border="0" cellpadding='0' cellspacing='0'>
       
       <tr>
         <td align="left">
           Version:
         </td>
         <td align="right">
           
           [&nbsp;<span class="var-sel">CMSSW_5_3_21</span>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_6_2_0_SLHC16">CMSSW_6_2_0_SLHC16</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_6_2_0_SLHC20">CMSSW_6_2_0_SLHC20</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_6_2_7">CMSSW_6_2_7</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_0_6">CMSSW_7_0_6</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_1_7">CMSSW_7_1_7</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_1_8">CMSSW_7_1_8</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_2_0_pre1">CMSSW_7_2_0_pre1</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_2_0_pre2">CMSSW_7_2_0_pre2</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_2_0_pre3">CMSSW_7_2_0_pre3</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_2_0_pre4">CMSSW_7_2_0_pre4</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_2_0_pre5">CMSSW_7_2_0_pre5</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_2_0_pre6">CMSSW_7_2_0_pre6</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_2_0_pre7">CMSSW_7_2_0_pre7</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_2_0_pre8">CMSSW_7_2_0_pre8</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py">CMSSW_7_2_1</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_3_0_pre1">CMSSW_7_3_0_pre1</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_3_X_2014-10-19-0200">CMSSW_7_3_X_2014-10-19-0200</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_3_X_2014-10-21-0200">CMSSW_7_3_X_2014-10-21-0200</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_3_X_2014-10-22-0200">CMSSW_7_3_X_2014-10-22-0200</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_3_X_2014-10-23-0200">CMSSW_7_3_X_2014-10-23-0200</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_3_X_2014-10-25-0200">CMSSW_7_3_X_2014-10-25-0200</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_3_X_2014-10-27-0200">CMSSW_7_3_X_2014-10-27-0200</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_3_X_2014-10-28-0200">CMSSW_7_3_X_2014-10-28-0200</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_3_X_2014-10-29-0200">CMSSW_7_3_X_2014-10-29-0200</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_3_X_2014-10-30-0200">CMSSW_7_3_X_2014-10-30-0200</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_2_X_2014-10-31-0200">CMSSW_7_2_X_2014-10-31-0200</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_3_X_2014-10-31-1400">CMSSW_7_3_X_2014-10-31-1400</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_3_X_2014-11-01-1400">CMSSW_7_3_X_2014-11-01-1400</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_1_1">CMSSW_7_1_1</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_3_X_2014-11-01-0200">CMSSW_7_3_X_2014-11-01-0200</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_2_X_2014-11-02-0200">CMSSW_7_2_X_2014-11-02-0200</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_3_X_2014-11-02-0200">CMSSW_7_3_X_2014-11-02-0200</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_1_X_2014-11-02-1400">CMSSW_7_1_X_2014-11-02-1400</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_3_X_2014-11-02-1400">CMSSW_7_3_X_2014-11-02-1400</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_3_X_2014-11-03-0200">CMSSW_7_3_X_2014-11-03-0200</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_3_X_2014-11-04-0200">CMSSW_7_3_X_2014-11-04-0200</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_2_X_2014-11-04-1400">CMSSW_7_2_X_2014-11-04-1400</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_3_X_2014-11-04-1400">CMSSW_7_3_X_2014-11-04-1400</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_3_X_2014-11-05-1400">CMSSW_7_3_X_2014-11-05-1400</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_3_X_2014-11-06-0200">CMSSW_7_3_X_2014-11-06-0200</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_1_0_pre10">CMSSW_7_1_0_pre10</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_2_1_patch2">CMSSW_7_2_1_patch2</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_3_X_2014-11-07-0200">CMSSW_7_3_X_2014-11-07-0200</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_2_X_2014-11-08-0200">CMSSW_7_2_X_2014-11-08-0200</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_3_X_2014-11-08-0200">CMSSW_7_3_X_2014-11-08-0200</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_1_X_2014-11-09-1400">CMSSW_7_1_X_2014-11-09-1400</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_2_X_2014-11-09-0200">CMSSW_7_2_X_2014-11-09-0200</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_3_X_2014-11-09-0200">CMSSW_7_3_X_2014-11-09-0200</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_3_X_2014-11-10-0200">CMSSW_7_3_X_2014-11-10-0200</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_3_X_2014-11-09-1400">CMSSW_7_3_X_2014-11-09-1400</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_3_X_2014-11-11-0200">CMSSW_7_3_X_2014-11-11-0200</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_3_0_pre2">CMSSW_7_3_0_pre2</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_2_1_patch3">CMSSW_7_2_1_patch3</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_2_1_patch4">CMSSW_7_2_1_patch4</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_3_X_2014-11-11-1400">CMSSW_7_3_X_2014-11-11-1400</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_2_X_2014-11-12-0200">CMSSW_7_2_X_2014-11-12-0200</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_3_X_2014-11-12-0200">CMSSW_7_3_X_2014-11-12-0200</a>&nbsp;]
           [&nbsp;<a class='varlink' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_7_3_X_2014-11-12-1400">CMSSW_7_3_X_2014-11-12-1400</a>&nbsp;]
           <br>
         </td>
       </tr>
     </table>
   </td>
  </tr>
</table>

<hr>
<pre class="file">
<a name=001 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#001">001</a> <a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=import">import</a> FWCore.<a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=ParameterSet">ParameterSet</a>.<a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=Config">Config</a> <a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=as">as</a> <a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=cms">cms</a>
<a name=002 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#002">002</a> 
<a name=003 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#003">003</a> <span class="comment"># configuration to model pileup for initial physics phase</span>
<a name=004 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#004">004</a> <span class="comment"></span><a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=from">from</a> SimGeneral.<a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=MixingModule">MixingModule</a>.mixObjects_cfi <a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=import">import</a> * 
<a name=005 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#005">005</a> <a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=from">from</a> SimGeneral.<a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=MixingModule">MixingModule</a>.mixPoolSource_cfi <a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=import">import</a> * 
<a name=006 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#006">006</a> 
<a name=007 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#007">007</a> <a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=mix">mix</a> = <a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=cms">cms</a>.<a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=EDProducer">EDProducer</a>(<span class='string'>"MixingModule"</span>,
<a name=008 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#008">008</a>     LabelPlayback = <a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=cms">cms</a>.<a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=string">string</a>(<span class='string'>''</span>),
<a name=009 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#009">009</a>     <a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=maxBunch">maxBunch</a> = <a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=cms">cms</a>.<a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=int32">int32</a>(3),
<a name=010 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#010">010</a>     <a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=minBunch">minBunch</a> = <a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=cms">cms</a>.<a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=int32">int32</a>(-12), <span class="comment">## in terms of 25 nsec</span>
<a name=011 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#011">011</a> <span class="comment"></span>
<a name=012 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#012">012</a>     bunchspace = <a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=cms">cms</a>.<a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=int32">int32</a>(50), <span class="comment">##ns</span>
<a name=013 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#013">013</a> <span class="comment"></span>    mixProdStep1 = <a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=cms">cms</a>.<a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=bool">bool</a>(<a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=False">False</a>),
<a name=014 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#014">014</a>     mixProdStep2 = <a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=cms">cms</a>.<a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=bool">bool</a>(<a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=False">False</a>),
<a name=015 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#015">015</a> 
<a name=016 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#016">016</a>     playback = <a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=cms">cms</a>.untracked.<a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=bool">bool</a>(<a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=False">False</a>),
<a name=017 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#017">017</a>     useCurrentProcessOnly = <a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=cms">cms</a>.<a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=bool">bool</a>(<a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=False">False</a>),
<a name=018 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#018">018</a>                    
<a name=019 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#019">019</a>     <a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=input">input</a> = <a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=cms">cms</a>.<a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=SecSource">SecSource</a>(<span class='string'>"PoolSource"</span>,
<a name=020 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#020">020</a>         <a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=type">type</a> = <a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=cms">cms</a>.<a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=string">string</a>(<span class='string'>'probFunction'</span>),
<a name=021 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#021">021</a>         nbPileupEvents = <a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=cms">cms</a>.<a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=PSet">PSet</a>(
<a name=022 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#022">022</a>           <a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=probFunctionVariable">probFunctionVariable</a> = <a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=cms">cms</a>.<a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=vint32">vint32</a>(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59),
<a name=023 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#023">023</a>           <a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=probValue">probValue</a> = <a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=cms">cms</a>.<a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=vdouble">vdouble</a>(
<a name=024 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#024">024</a>                   2.90E-15,
<a name=025 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#025">025</a>                   2.00E-09,
<a name=026 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#026">026</a>                   1.57E-06,
<a name=027 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#027">027</a>                   1.46E-04,
<a name=028 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#028">028</a>                   3.10E-04,
<a name=029 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#029">029</a>                   4.83E-05,
<a name=030 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#030">030</a>                   5.74E-05,
<a name=031 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#031">031</a>                   5.36E-05,
<a name=032 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#032">032</a>                   4.45E-04,
<a name=033 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#033">033</a>                   4.33E-03,
<a name=034 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#034">034</a>                   2.00E-02,
<a name=035 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#035">035</a>                   4.26E-02,
<a name=036 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#036">036</a>                   6.00E-02,
<a name=037 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#037">037</a>                   7.35E-02,
<a name=038 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#038">038</a>                   8.63E-02,
<a name=039 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#039">039</a>                   9.54E-02,
<a name=040 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#040">040</a>                   9.82E-02,
<a name=041 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#041">041</a>                   9.41E-02,
<a name=042 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#042">042</a>                   8.57E-02,
<a name=043 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#043">043</a>                   7.62E-02,
<a name=044 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#044">044</a>                   6.66E-02,
<a name=045 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#045">045</a>                   5.64E-02,
<a name=046 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#046">046</a>                   4.57E-02,
<a name=047 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#047">047</a>                   3.51E-02,
<a name=048 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#048">048</a>                   2.49E-02,
<a name=049 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#049">049</a>                   1.59E-02,
<a name=050 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#050">050</a>                   9.08E-03,
<a name=051 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#051">051</a>                   4.74E-03,
<a name=052 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#052">052</a>                   2.30E-03,
<a name=053 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#053">053</a>                   1.07E-03,
<a name=054 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#054">054</a>                   4.78E-04,
<a name=055 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#055">055</a>                   2.06E-04,
<a name=056 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#056">056</a>                   8.52E-05,
<a name=057 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#057">057</a>                   3.32E-05,
<a name=058 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#058">058</a>                   1.20E-05,
<a name=059 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#059">059</a>                   3.99E-06,
<a name=060 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#060">060</a>                   1.21E-06,
<a name=061 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#061">061</a>                   3.29E-07,
<a name=062 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#062">062</a>                   8.05E-08,
<a name=063 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#063">063</a>                   1.76E-08,
<a name=064 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#064">064</a>                   3.45E-09,
<a name=065 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#065">065</a>                   5.99E-10,
<a name=066 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#066">066</a>                   9.26E-11,
<a name=067 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#067">067</a>                   1.27E-11,
<a name=068 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#068">068</a>                   1.54E-12,
<a name=069 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#069">069</a>                   1.66E-13,
<a name=070 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#070">070</a>                   1.58E-14,
<a name=071 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#071">071</a>                   1.33E-15,
<a name=072 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#072">072</a>                   9.91E-17,
<a name=073 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#073">073</a>                   6.53E-18,
<a name=074 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#074">074</a>                   3.84E-19,
<a name=075 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#075">075</a>                   1.65E-20,
<a name=076 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#076">076</a>                   0.00E+00,
<a name=077 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#077">077</a>                   0.00E+00,
<a name=078 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#078">078</a>                   0.00E+00,
<a name=079 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#079">079</a>                   0.00E+00,
<a name=080 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#080">080</a>                   0.00E+00,
<a name=081 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#081">081</a>                   0.00E+00,
<a name=082 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#082">082</a>                   0.00E+00,
<a name=083 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#083">083</a>                   0.00E+00
<a name=084 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#084">084</a>                 ),
<a name=085 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#085">085</a>           <a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=histoFileName">histoFileName</a> = <a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=cms">cms</a>.untracked.<a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=string">string</a>(<span class='string'>'histProbFunction.root'</span>),
<a name=086 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#086">086</a>         ),
<a name=087 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#087">087</a>     sequential = <a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=cms">cms</a>.untracked.<a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=bool">bool</a>(<a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=False">False</a>),                          
<a name=088 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#088">088</a>         manage_OOT = <a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=cms">cms</a>.untracked.<a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=bool">bool</a>(<a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=True">True</a>),  <span class="comment">## manage out-of-time pileup</span>
<a name=089 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#089">089</a> <span class="comment"></span>        <span class="comment">## setting this to True means that the out-of-time pileup</span>
<a name=090 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#090">090</a> <span class="comment"></span>        <span class="comment">## will have a different distribution than in-time, given</span>
<a name=091 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#091">091</a> <span class="comment"></span>        <span class="comment">## by what is described on the next line:</span>
<a name=092 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#092">092</a> <span class="comment"></span>        OOT_type = <a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=cms">cms</a>.untracked.<a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=string">string</a>(<span class='string'>'Poisson'</span>),  <span class="comment">## generate OOT with a Poisson matching the number chosen for in-time</span>
<a name=093 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#093">093</a> <span class="comment"></span>        <span class="comment">#OOT_type = cms.untracked.string('fixed'),  ## generate OOT with a fixed distribution</span>
<a name=094 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#094">094</a> <span class="comment"></span>        <span class="comment">#intFixed_OOT = cms.untracked.int32(2),</span>
<a name=095 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#095">095</a> <span class="comment"></span>        <a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=fileNames">fileNames</a> = FileNames 
<a name=096 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#096">096</a>     ),
<a name=097 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#097">097</a>     mixObjects = <a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=cms">cms</a>.<a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=PSet">PSet</a>(
<a name=098 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#098">098</a>         mixCH = <a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=cms">cms</a>.<a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=PSet">PSet</a>(
<a name=099 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#099">099</a>             mixCaloHits
<a name=100 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#100">100</a>         ),
<a name=101 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#101">101</a>         mixTracks = <a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=cms">cms</a>.<a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=PSet">PSet</a>(
<a name=102 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#102">102</a>             mixSimTracks
<a name=103 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#103">103</a>         ),
<a name=104 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#104">104</a>         mixVertices = <a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=cms">cms</a>.<a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=PSet">PSet</a>(
<a name=105 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#105">105</a>             mixSimVertices
<a name=106 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#106">106</a>         ),
<a name=107 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#107">107</a>         mixSH = <a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=cms">cms</a>.<a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=PSet">PSet</a>(
<a name=108 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#108">108</a>             mixSimHits
<a name=109 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#109">109</a>         ),
<a name=110 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#110">110</a>         mixHepMC = <a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=cms">cms</a>.<a class='fid' href="/SDT/lxr/ident?v=CMSSW_5_3_21;i=PSet">PSet</a>(
<a name=111 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#111">111</a>             mixHepMCProducts
<a name=112 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#112">112</a>         )
<a name=113 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#113">113</a>     )
<a name=114 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#114">114</a> )
<a name=115 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#115">115</a> 
<a name=116 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#116">116</a> 
<a name=117 class='fline' href="/SDT/lxr/source/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21#117">117</a> </pre>
<hr>
<table width="100%" cellpadding="0" border="0">
  <tr valign="middle">
    
    <td align="center" nowrap="nowrap">
      [&nbsp;<span class='modes-sel'>source navigation</span>&nbsp;]</td>
    <td align="center" nowrap="nowrap">
      [&nbsp;<a class='modes' href="/SDT/lxr/diff/SimGeneral/MixingModule/python/mix_2012A_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21">diff markup</a>&nbsp;]</td>
    <td align="center" nowrap="nowrap">
      [&nbsp;<a class='modes' href="/SDT/lxr/ident?v=CMSSW_5_3_21">identifier search</a>&nbsp;]</td>
    <td align="center" nowrap="nowrap">
      [&nbsp;<a class="modes" href="/SDT/lxr/search?v=CMSSW_5_3_21">general search</a>&nbsp;]</td>
  </tr>
</table>
<hr>
<table width="100%" cellpadding="0" border="0">
  <tr>
    <td align="left">
      This page was automatically generated by the 
      <a href="http://lxr.sf.net/">LXR engine</a>.
      <address>
        <a href="mailto:lxr-general@lists.sf.net">The LXR team</a>
      </address>
    </td>
    <td align="right">
      <a href="http://validator.w3.org/check/referer"><img border="0"
        src="templates/valid-html401.png"
        alt="Valid HTML 4.01!" height="31" width="88"></a>
    </td>
  </tr>
</table>

</body>
</html>
