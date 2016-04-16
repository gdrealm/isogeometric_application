//problemtype for Isogeometric simulation module for Kratos
//(c) 2013 Hoang Giang Bui, Ruhr-University Bochum


*set var i=0
*set var j=tcl(GetNumPoints)
*for(i=1;i<=j;i=i+1) *\
*set var num=i
*tcl(WritePointInfoVerbose *num)

*end


*set var i=0
*set var j=tcl(GetNumLines)
*for(i=1;i<=j;i=i+1) *\
*set var num=i
*tcl(WriteLineInfoVerbose *num)

*end


*set var i=0
*set var j=tcl(GetNumSurfaces)
*for(i=1;i<=j;i=i+1) *\
*set var num=i
*tcl(WriteSurfaceInfoVerbose *num)

*end


*set var i=0
*set var j=tcl(GetNumVolumes)
*for(i=1;i<=j;i=i+1) *\
*set var num=i
*tcl(WriteVolumeInfoVerbose *num)

*end


