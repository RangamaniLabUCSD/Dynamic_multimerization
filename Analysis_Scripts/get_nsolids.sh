while getopts f: flag
do
    case "${flag}" in
        f) filepath=${OPTARG};;
    esac
done
filename="${filepath}/objects.cmo"
startl=$(grep -n -a "section solid" $filename|head -1 | cut -d : -f 1)
endlA=$(grep -n -a "section single F" $filename|head -1 | cut -d : -f 1)
endlB=$(grep -n -a "section single A" $filename|head -1 | cut -d : -f 1)
endl=$endlA
if [ "$endlA" -gt "$endlB" ] && [ "$endlB" -gt 0]
   then endl=$endlB
fi
var=5 
nsolids="$(((endl-startl-1)/var))"
echo $nsolids