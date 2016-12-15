BEGIN {
FS="@";
}

{
printf("%s (\"%s\", %s, %s),\n", $2, $3, $1,$4) 
next
}