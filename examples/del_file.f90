open(unit=1234, iostat=stat, file=file, status='old')
if (stat == 0) close(1234, status='delete')