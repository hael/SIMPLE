#ifndef THROW_HARD
#define THROW_HARD(msg) call simple_exception(msg, __FILENAME__ , __LINE__)
#endif
#ifndef THROW_WARN
#define THROW_WARN(msg) call simple_exception(msg, __FILENAME__ , __LINE__, l_stop=.false.)
#endif
