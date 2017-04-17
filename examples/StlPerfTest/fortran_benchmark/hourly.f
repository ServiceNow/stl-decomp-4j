c data are 89345 hourly values (10 years) in file hourly_stl_test.csv
      real y(89345), season(89345), trend(89345), rw(89345)
      real work(106817,7)
      real start, finish
      logical robust
      read(5,*)(y(i), i=1, 89345)
      n=89345
      np=8736
      ns=893451
      nt=13105
      nl=8737
      no=0
      ni=2
      nsjump=89346
      ntjump=1311
      nljump=874
      isdeg=0
      itdeg=1
      ildeg=1
      iterations = 100

      call cpu_time(start)

      do 1001 i = 1, iterations

         call stl(y,n,np,ns,nt,nl,isdeg,itdeg,ildeg,nsjump,ntjump,
     &        nljump,ni,no,rw,season, trend, work)

 1001 continue

      call cpu_time(finish)

      print '("Time = ",f6.3," ms.")', 1000.0 * (finish-start)
     &     / iterations

      open(unit=99, file="output.csv", status='replace')

      do 2002 i = 1, 89345
         d = y(i)
         s = season(i)
         t = trend(i)
         r = y(i) - season(i) - trend(i)
         write(unit=99,  FMT="(1x, E24.16, (3(',',E24.16)))") d, s, t, r
 2002 continue

      close(unit=99)
      
      end
