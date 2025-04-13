#include <pthread.h>
// https://docs.oracle.com/cd/E19455-01/806-5257/6je9h033p/index.html#:~:text=Solaris%20threads%20supports%20functions%20that,content%20is%20effectively%20the%20same.
#define THR_BOUND PTHREAD_SCOPE_SYSTEM
#define THR_DAEMON 0
#define thr_create(base, size, wrkr, arg, attr, tid)  pthread_create(tid, attr, wrkr, arg)
#define thr_yield   sched_yield // POSIX.4
#if 0
#define thr_exit    pthread_exit
#define thr_join    pthread_join
#define thr_self    pthread_self
#define thr_kill    pthread_kill
#define thr_sigsetmask      pthread_sigmask
#define thr_setprio         pthread_setschedparam
#define thr_getprio         pthread_getschedparam
#define thr_setconcurrency  pthread_setconcurrency
#define thr_getconcurrency  pthread_getconcurrency
#define thr_suspend   - 
#define thr_continue  - 
#define thr_keycreate       pthread_key_create
-                     pthread_key_delete
#define thr_setspecific     pthread_setspecific
#define thr_getspecific     pthread_getspecific
-                     pthread_once
-                     pthread_equal
-                     pthread_cancel
-                     pthread_testcancel
-                     pthread_cleanup_push
-                     pthread_cleanup_pop
-                     pthread_setcanceltype
-                     pthread_setcancelstate
#endif
//
//
typedef pthread_mutex_t mutex_t;
#define USYNC_THREAD  NULL // PTHREAD_PROCESS_PRIVATE
#define mutex_init(a,b,c)    pthread_mutex_init((a), (b))
#define mutex_lock    pthread_mutex_lock
#define mutex_unlock  pthread_mutex_unlock
#if 0
#define mutex_trylock  pthread_mutex_trylock
#define mutex_destroy  pthread_mutex_destroy
#endif
//
typedef pthread_cond_t  cond_t;
#define cond_init(a, b, c)       pthread_cond_init((a), (b))
#define cond_wait       pthread_cond_wait
#define cond_broadcast  pthread_cond_broadcast
#define cond_signal     pthread_cond_signal
#if 0
#define cond_timedwait  pthread_cond_timedwait
#define cond_destroy    pthread_cond_destroy
#endif
//
#if 0
rwlock_init     pthread_rwlock_init
rwlock_destroy  pthread_rwlock_destroy
rw_rdlock       pthread_rwlock_rdlock
rw_wrlock       pthread_rwlock_wrlock
rw_unlock       pthread_rwlock_unlock
rw_tryrdlock    pthread_rwlock_tryrdlock
rw_trywrlock    pthread_rwlock_trywrlock
-   pthread_rwlockattr_init
-   pthread_rwlockattr_destroy
-   pthread_rwlockattr_getpshared
-   pthread_rwlockattr_setpshared

// semaphores
sema_init     sem_init POSIX 1003.4
sema_destroy  sem_destroy POSIX 1003.4
sema_wait     sem_wait POSIX 1003.4
sema_post     sem_post POSIX 1003.4

sema_trywait  sem_trywait POSIX 1003.4

// fork
fork1         fork
-               pthread_atfork

fork (multiple thread copy) - 

-   pthread_mutexattr_init
-   pthread_mutexattr_destroy

type argument in cond_init  pthread_mutexattr_setpshared

-   pthread_mutexattr_getpshared

-   pthread_mutex_attr_settype

-   pthread_mutex_attr_gettype

-   pthread_condattr_init

-   pthread_condattr_destroy

type argument in cond_init  pthread_condattr_setpshared

-   pthread_condattr_getpshared

-   pthread_attr_init

-   pthread_attr_destroy

THR_BOUND flag in thr_create  pthread_attr_setscope

-   pthread_attr_getscope

-   pthread_attr_setguardsize

-   pthread_attr_getguardsize

stack_size argument in thr_create  pthread_attr_setstacksize

-   pthread_attr_getstacksize

stack_addr argument in thr_create  pthread_attr_setstackaddr

-   pthread_attr_getstackaddr

THR_DETACH flag in thr_create  pthread_attr_setdetachstate

-   pthread_attr_getdetachstate

-   pthread_attr_setschedparam

-   pthread_attr_getschedparam

-   pthread_attr_setinheritsched

-   pthread_attr_getinheritsched

-   pthread_attr_setsschedpolicy

-   pthread_attr_getschedpolicy
#endif
