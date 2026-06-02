# Useful GDB Commands

```shell
gdb -ex=r –args myprogram arg1 arg2
```
`-ex=r` is short for `-ex=run` and tells gdb to run your program immediately, rather than wait for you to type "run" at the prompt. Then `--args` says that everything that follows is the command and arguments, just as you'd normally type them at the commandline prompt.