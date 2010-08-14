/*
 *  This file is part of Distributed Replica.
 *  Copyright May 9 2009
 *
 *  Distributed Replica manages a series of simulations that separately sample phase space
 *  and coordinates their efforts under the Distributed Replica Potential Energy Function.
 *  See, for example T. Rodinger, P.L. Howell, and R. Pomès, "Distributed Replica Sampling"    
 *  J. Chem. Theory Comput., 2:725 (2006).
 *
 *  Distributed Replica is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Distributed Replica is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Distributed Replica.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <string>

using namespace std;

double string_double(string s);

double string_double(string s)
{
 int p=0,ps=0;
 double num=0.0;
 int sign=1;
 bool dot=true,sc=true;

 if(s[0]=='-') //negative number
   sign=-1;

 for(unsigned int i=0; i<s.size(); ++i){ // find e or E
   if(s[i]=='e'||s[i]=='E'){ ps=i; break; }
   if(i==s.size()-1) {ps=i+1; sc=false; }
 }

 for(unsigned int i=0; i<s.size(); ++i){ // find .
   if(s[i]=='.'){ p=i; break; }
   if(i==s.size()-1) {p=i+1; dot=false; }
 }

 if(sc){ 
   double c=1.0;
   if(!dot){
     for(int i=ps-1; i>0; --i){
       num+=c*(s[i]-48);
       c*=10.0;
     }
     if(sign==1) num+=c*(s[0]-48);
   }
   else{
     for(int i=p-1; i>0; --i){
       num+=c*(s[i]-48);
       c*=10.0;
     }
     if(sign==1) num+=c*(s[0]-48);
     c=0.1;
     for(int i=p+1; i<ps; ++i){
       num+=c*(s[i]-48);
       c*=0.1;
     }
   }
   int exp=0,ci=1;
   for(int i=s.size()-1; i>ps+1; --i){
     exp+=ci*(s[i]-48);
     ci*=10;
   }
   if(s[ps+1]=='-'){
     for(int i=0; i<exp; ++i)
       num*=0.1;
   }
   if(s[ps+1]=='+'){
     for(int i=0; i<exp; ++i)
       num*=10.0;
   }
   else{
     exp+=ci*(s[ps+1]-48);
     for(int i=0; i<exp; ++i)
       num*=10.0;
   }
 }
 else{
   double c=1.0;
   for(int i=p-1; i>0; --i){
     num+=c*(s[i]-48);
     c*=10.0;
   }
   if(sign==1) num+=c*(s[0]-48);

   if(dot){
     c=0.1;
     for(unsigned int i=p+1; i<s.size(); ++i){
       num+=c*(s[i]-48);
       c*=0.1;
     }
   }
 }
 return sign*num;
}
