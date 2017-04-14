c data are 708 monthly values (59 years) in file co2_2017
      real y(708), season(708), trend(708), rw(708), work(732,7)
      real start, finish
      logical robust
      read(5,*)(y(i), i=1, 708)
      n=708
      np=12
      ns=35
      nt=19
      nl=13
      no=0
      ni=2
      nsjump=4
      ntjump=2
      nljump=2
      isdeg=1
      itdeg=1
      ildeg=1
      iterations = 10000

      call cpu_time(start)

      do 1001 i = 1, iterations

         call stl(y,n,np,ns,nt,nl,isdeg,itdeg,ildeg,nsjump,ntjump,
     &        nljump,ni,no,rw,season, trend, work)

 1001 continue

      call cpu_time(finish)

      print '("Time = ",f6.3," ms.")', 1000.0 * (finish-start)
     &     / iterations

      open(unit=99, file="output.csv", status='replace')

      do 2002 i = 1, 708
         d = y(i)
         s = season(i)
         t = trend(i)
         r = y(i) - season(i) - trend(i)
         write(unit=99,  FMT="(1x, E24.16, (3(',',E24.16)))") d, s, t, r
 2002 continue

      close(unit=99)
      
      end
