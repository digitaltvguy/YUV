##YUV


Conversion between RGB and YUV formats. Supports: 709, 2020, XYZ.  Additionally supports 2020 constant luminance (for just 1k and 10k PQ nits)

Requirements: libtiff

Usage: 
tif2yuv  `<frame name>`  `<resolution>`  `<bit depth>`  `<transfer function>`


yuv2tiff `<yuv file name>` `<resolution>` `<bit depth>` `<transfer function>` `<# frames>`

`<frame name>` must be tiff file written as one row per strip 16 bits RGB order. (ctlrender does this)

`<resolution>` default is 3840x2160 other options are HD1920 (for 1920x1080).  It will try to center cut if input is wide but untested. Safest is to use exactly 3840x2160 or 1920x1080 rasters only.

`<bit depth>` 12 bits is default does not need to be signalled. Other options are B10 and B14 for 10 and 12 bits. Acheived by right shifting into the final 16 bits written into yuv and left shifting reading out from yuv back to tiff.

`<transfer function>` options are "709" or "2020" and work as expected. XYZ is the default and will produce a YUV file in YDzDx format. "Y500" is the compensated YDzDx processing that is in the SMPTE 2084 specification. It is only reccommended to use this setting.
"2020CL" is for constant luminance this option must also know Gamma and only "G1k" (1k nit PQ) or "G10k" (10k nit PQ) is currently supported.

`<# frames>` `"-f #"`  only used with yuv2tiff to tell it to only write out that many frames.

Example:

`# writing a 2020 Constant Luminance HD yuv file:`

$EDRHOME/Tools/YUV/tif2yuv $framename B10 2020C G1k HD1920


`# writing out only the first 5 frames of YDzDx SMPTE 2084`

`# 12 bits is default, 3840 is default`

$EDRHOME/Tools/YUV/tif2yuv $yuvfile Y500 -f 5



Reccommended directory structure:

$EDRHOME : base directory for all projects, contains folders for utilities and further testing

$EDRDATA = $EDRHOME/EDRDATA : base directory for testing

$EDRHOME/ACES : base directory for AMPAS ctlrender builds, IlmLib, OpenEXR and CTL scripts

$EDRHOME/DCP : base directory for DCP testing

$EDRHOME/FF : base directory for FFMPEG builds (if needed).

$EDRHOME/HEVC : base directory for HM and x265

$EDRHOME/Tools : base directory for tools and other code (YUV, pattern, tifcmp, sigma_compare etc..)



