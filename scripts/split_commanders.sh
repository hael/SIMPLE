#!/bin/bash

infile=$1 #'../src/main/simple_commander_misc.f90'
base=$(basename "${infile}" .f90)

rm -f tmp.txt template.header ${base}* local.tmp sub.tmp type.tmp external.tmp
sed -n '1h;1!H;${;g;s/subroutine \([_a-z0-9A-Z]*\)(\([^(]*\))\(.*\)end subroutine \1\s/\nXXX\1\n subroutine \1(\2)\n\3\nend subroutine \1\nYYY\n/g;p;}' "${infile}" > tmp.txt

# sed -n '1h;1!H;${;g;s/subroutine \([_a-z]*\)(\([^(]*\))\(.*\)end subroutine \1/XXX\nname:\1\nargs:\2 \nsub:\3 YYY\n/g;p;}' ../src/main/simple_commander_misc.f90

sed -n '1h;1!H;${;g;s/module \(.*\)include\(.*\)implicit none.*/module \1\n include \2\nimplicit none/g;p;}' "${infile}" | tail -n +3 > template.header


for i in $(grep XXX tmp.txt| sed -e 's/XXXexec_\(.*\)$/\1/' ); do
    rm -f "${base}_${i}.f90" type.tmp sub.tmp local.tmp external.tmp
    grep -B2 -A1 -e "=>[ ]*exec_"${i}"\b" < tmp.txt > type.tmp
    sed -n '/XXXexec_'"${i}"'\b/{:start /YYY/!{N;b start};s/XXXexec_'"${i}"'\b\(.*\)YYY/\1/p}' < tmp.txt > sub.tmp
    grep -E  '^[ ]*type(.*_commander)' sub.tmp  > local.tmp
    if [ -s local.tmp ];then
        while read -r internalcommander ; do
            shortcommander=$(echo "${internalcommander}" | sed 's/^[ ]*type(\(.*\)_commander).*/\1/');
            /bin/grep -E -q  "^[ ]*use simple.*only:.*${shortcommander}_commander" sub.tmp
            [ $? == 1 ] && echo "${shortcommander}"
        done < local.tmp > external.tmp
        sed -i 's/^\([^!].*\)/use \'"${base}"'_\1, only: \1_commander/'  external.tmp
        [ -s external.tmp ] && sed -i '3r external.tmp' sub.tmp
    fi

    (echo "module ${base}_${i}"
     cat template.header
     echo public :: ${i}_commander
     echo "private"
     echo "#include \"simple_local_flags.inc\""
     cat type.tmp
     echo "contains"
     cat sub.tmp
     echo ""
     echo "end module "
    ) > "${base}_${i}.f90"
done

(
echo "module $base"
for i in $(grep XXX tmp.txt| sed -e 's/XXXexec_\(.*\)$/\1/' ); do
    echo "use  ${base}_${i}"
done
echo "implicit none"
echo  "contains"
echo "end module"
)> "$base.f90"

rm -f type.tmp sub.tmp template.header tmp.txt local.tmp external.tmp
