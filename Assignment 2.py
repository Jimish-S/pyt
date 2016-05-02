class HeatExchanger:
      mA=1.0 #mass flow rate of fluid A
      mB=2.0 #mass flow rate of fluid B
      U=100.0 #overall heat transfer coefficient
      P=1.0 #Perimeter of heat exchanger over which heat tranfer takes place
      L=10.0 #length of the heat exchanger
      Tain=400 #inlet temperature of fluid A
      Tbin=300 #inlet temperature of fluid B
      n=10 #no of discretizations
          
      def Heat(self):
          mA=self.mA
          mB=self.mB
          U=self.U
          P=self.P
          L=self.L
          Tain=self.Tain
          Tbin=self.Tbin
          n=self.n
          deltal=L/(n-1)
          n_list=range(n)
          j=1
          e=1000;
          while e>0.001:
              Ta1=Tain
              Tb1=Tbin+j*0.01
              for i in n_list:
                  Ta= Ta1
                  Tb= Tb1
                  def Cpa(T):
                      Cpa= 4000.0 + 10.0*T+ 0.01*T*T
                      return Cpa
                  def Cpb(T):
                      Cpb= 3000.0 +5.0*T + 0.02*T*T
                      return Cpb
                  ca=Cpa(Ta)
                  cb=Cpb(Tb)
                  dtadl=-((U*P)*(Ta-Tb)/(mA*ca))
                  dtbdl=-((U*P)*(Ta-Tb)/(mB*cb))
                  Ta1=(dtadl*deltal)+Ta
                  Tb1=(dtbdl*deltal)+Tb
              e=Tb-Tbin   
              if e<0:
                  e=Tbin-Tb
              j=j+1;    
          print "outlet temperature of B", Tb;
          print "outlet temperature of A", Ta;
          return Tb          
                  
                  
                  
              
             
              
              

