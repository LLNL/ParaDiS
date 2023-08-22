
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include "String_t.h"

//---------------------------------------------------------------------------------------------------------

String_t::String_t(void) { str=0; }

String_t::String_t(const char *s0, const char *s1, const char *s2, const char *s3, const char *s4,
                   const char *s5, const char *s6, const char *s7, const char *s8, const char *s9 )
{
   str=0; Append(s0,s1,s2,s3,s4,s5,s6,s7,s8,s9);
}

String_t::String_t(const String_t & s)    { str = ( s.str ? strdup(s.str) : 0 ); }

//---------------------------------------------------------------------------------------------------------

String_t::~String_t() { if (str) delete [] str; str=0; }

//---------------------------------------------------------------------------------------------------------

const String_t & String_t::operator = (const String_t & s) 
{
   if (this!=&s)
   {
      if (str) delete [] str;  
      str = ( s.str ? strdup(s.str) : 0 ); 
   }

   return(*this);
}

//---------------------------------------------------------------------------------------------------------

const String_t & String_t::operator = (const char *s) 
{ 
   if (str) delete [] str;  
   str = ( s ? strdup(s) : 0 ); 

   return(*this); 
}

//---------------------------------------------------------------------------------------------------------

int String_t::Length() const { return( str ? strlen(str) : 0 ); }

static long long str_ll (const char *str)
{
   long long s = 0;

   char *p  = (char *) str;
   char *pe = (char *) str + strlen(str);

   if (p && pe && (p<pe))
   {
      if (isspace(*p)) { while( (p<pe) && (isspace(*p)) ) p++; }    // skip leading whitespace

      int base = 10; 

      if (p<pe)
      {
         // if the leading digit is zero, the string contains a hex, binary, or octal number...

         if (p[0]=='0')
         {
                 if ( (p[1]=='x') || (p[1]=='X') ) { base = 16; }   // (hex)
            else if ( (p[1]=='b') || (p[1]=='B') ) { base =  2; }   // (binary)
            else                                   { base =  8; }   // (octal)
         }
      }

      if (base==2)
      {
         p++; p++;
         while( (p<pe) && ((*p=='0') || (*p=='1')) ) { s = (s<<1) | ( (*p=='1') ? 1 : 0 ); p++; }
      }
      else
         s = strtoll(p,0,base);
   }

   return(s);
}

//---------------------------------------------------------------------------------------------------------

String_t::operator          char     *() const { return( str ); }
String_t::operator          int       () const { return( str ? (         int      ) str_ll(str)   : 0    ); }
String_t::operator unsigned int       () const { return( str ? (unsigned int      ) str_ll(str)   : 0    ); }
String_t::operator          short     () const { return( str ? (         short    ) str_ll(str)   : 0    ); }
String_t::operator unsigned short     () const { return( str ? (unsigned short    ) str_ll(str)   : 0    ); }
String_t::operator          long      () const { return( str ? (         long     ) str_ll(str)   : 0    ); }
String_t::operator unsigned long      () const { return( str ? (unsigned long     ) str_ll(str)   : 0    ); }
String_t::operator          long long () const { return( str ? (         long long) str_ll(str)   : 0    ); }
String_t::operator unsigned long long () const { return( str ? (unsigned long long) str_ll(str)   : 0    ); }
String_t::operator          float     () const { return( str ? (         float    ) strtod(str,0) : 0.0f ); }
String_t::operator          double    () const { return( str ? (         double   ) strtod(str,0) : 0.0  ); }

//---------------------------------------------------------------------------------------------------------

void String_t::operator += (const char *s) 
{
   if (s)
   {
      int n  = ( str ? strlen(str) : 0 );
          n += ( s   ? strlen(s  ) : 0 );

      char *tmp = ( (n>0) ? new char[n+1] : 0 );

      if (tmp)
      {
         tmp[0]=0;
 
         if (str) strcat(tmp,str);
         if (s  ) strcat(tmp,s  );
      }

      if (str) delete [] str;

      str = tmp;
   }
}

void String_t::operator += (const String_t  &s) { *this += (s.str); }

//---------------------------------------------------------------------------------------------------------

int String_t::operator == (const char    *s) const { return( (str && s    ) ? ( (strcmp(str,s    )==0) ? 1 : 0 ) : 0 ); }
int String_t::operator != (const char    *s) const { return( (str && s    ) ? ( (strcmp(str,s    )!=0) ? 1 : 0 ) : 0 ); }
int String_t::operator <  (const char    *s) const { return( (str && s    ) ? ( (strcmp(str,s    )< 0) ? 1 : 0 ) : 0 ); }
int String_t::operator >  (const char    *s) const { return( (str && s    ) ? ( (strcmp(str,s    )> 0) ? 1 : 0 ) : 0 ); }
int String_t::operator == (const String_t  &s) const { return( (str && s.str) ? ( (strcmp(str,s.str)==0) ? 1 : 0 ) : 0 ); }
int String_t::operator != (const String_t  &s) const { return( (str && s.str) ? ( (strcmp(str,s.str)!=0) ? 1 : 0 ) : 0 ); }
int String_t::operator <  (const String_t  &s) const { return( (str && s.str) ? ( (strcmp(str,s.str)< 0) ? 1 : 0 ) : 0 ); }
int String_t::operator >  (const String_t  &s) const { return( (str && s.str) ? ( (strcmp(str,s.str)> 0) ? 1 : 0 ) : 0 ); }

//---------------------------------------------------------------------------------------------------------

int  String_t::Equiv (const char *s) const
{
   char *p   = ( str ? (char *) str : 0 );
   char *ps  = ( s   ? (char *) s   : 0 );

   if (p && ps)
   {
      if ( strlen(p) != strlen(ps) ) return(0);

      char *pe = p+strlen(p);

      while (p<pe) { if ( tolower(*p++) != tolower(*ps++) ) return(0); }
   }

   return(0);
}

//---------------------------------------------------------------------------------------------------------

void String_t::Append 
(
   const char *s0, const char *s1, const char *s2, const char *s3, const char *s4,
   const char *s5, const char *s6, const char *s7, const char *s8, const char *s9 
)
{
   int n =0;

   n += ( str ? strlen(str) : 0 );
   n += ( s0  ? strlen(s0)  : 0 );
   n += ( s1  ? strlen(s1)  : 0 );
   n += ( s2  ? strlen(s2)  : 0 );
   n += ( s3  ? strlen(s3)  : 0 );
   n += ( s4  ? strlen(s4)  : 0 );
   n += ( s5  ? strlen(s5)  : 0 );
   n += ( s6  ? strlen(s6)  : 0 );
   n += ( s7  ? strlen(s7)  : 0 );
   n += ( s8  ? strlen(s8)  : 0 );
   n += ( s9  ? strlen(s9)  : 0 );

   char *tmp = ( (n>0) ? new char[n+1] : 0 );

   if (tmp)
   {
      memset(tmp,0,(n+1)*sizeof(char));

      if (str) strcat(tmp,str);
      if (s0 ) strcat(tmp,s0 );
      if (s1 ) strcat(tmp,s1 );
      if (s2 ) strcat(tmp,s2 );
      if (s3 ) strcat(tmp,s3 );
      if (s4 ) strcat(tmp,s4 );
      if (s5 ) strcat(tmp,s5 );
      if (s6 ) strcat(tmp,s6 );
      if (s7 ) strcat(tmp,s7 );
      if (s8 ) strcat(tmp,s8 );
      if (s9 ) strcat(tmp,s9 );
   }

   if (str) { delete [] str; str=0; }

   str = tmp;
}

//---------------------------------------------------------------------------------------------------------

void String_t::Upper (void) { if (str) { int n=(int) strlen(str); for (int i=0; (i<n); ++i) str[i]=toupper(str[i]); } }
void String_t::Lower (void) { if (str) { int n=(int) strlen(str); for (int i=0; (i<n); ++i) str[i]=tolower(str[i]); } }

//---------------------------------------------------------------------------------------------------------

void String_t::Replace (const char c0, const char c1) 
{ 
   if (str) 
   { 
      for (int i=0, n=strlen(str); (i<n); ++i) 
         str[i] = ( str[i]==c1 ? c0 : str[i] );
   } 
}

//---------------------------------------------------------------------------------------------------------

void String_t::Strip (void)
{
   Strip_L();
   Strip_R();
}

//---------------------------------------------------------------------------------------------------------

void String_t::Strip_L (void)
{
   if (str)
   {
      char *p  = str;
      char *pe = str + strlen(str);
      char *q  = str;

      while ( (p<pe) && (isspace(*p)) ) { p++; }

      if (q<p)
      {
         while ( (p<pe) ) { *q++ = *p++; }
         while ( (q<pe) ) { *q++ = 0;    }
      }
   }
}

//---------------------------------------------------------------------------------------------------------

void String_t::Strip_R (void)
{
   if (str)
   {
      char *p  = str;
      char *pe = str + strlen(str);

      while ( (p<pe) && (isspace(*(pe-1))) ) { *(pe-1)=0; pe--; } 
   }
}

//---------------------------------------------------------------------------------------------------------

int String_t::Prefix (const char *pfx) const
{
   if (pfx && str)
   {
      int n = strlen(pfx);
      int m = strlen(str);

      if ( (n>0) && (n<=m) )
         return ( (strncmp(str,pfx,n)==0) ? 1 : 0 );
   }

   return(0);
}

//---------------------------------------------------------------------------------------------------------

int String_t::Suffix (const char *sfx) const
{
   if (sfx && str)
   {
      int n = strlen(sfx);
      int m = strlen(str);

      if ( (n>0) && (n<=m) )
         return ( (strncmp(str+m-n,sfx,n)==0) ? 1 : 0 );
   }

   return(0);
}

//---------------------------------------------------------------------------------------------------------

char *String_t::Find (const char c) const
{
   if (str)
   {
      char *p  = str;
      char *pe = str+strlen(str);

      while ( (p<pe) && (*p!=c) ) { p++; }

      return ( (*p==c) ? p : 0 );
   }

   return(0);
}

//---------------------------------------------------------------------------------------------------------

char *String_t::Find (const char *s) const
{
   if (str && s)
   {
      int   n  = strlen(s);
      char *p  = (char *) str ;
      char *pe = (char *) str+strlen(str)-n;
   
      if (p && pe && (p<=pe))
      {      
         while ((p<=pe) && (strncmp(p,s,n)!=0)) { p++; }
  
         return ( (p && (p<=pe) && (strncmp(p,s,n)==0)) ? p : 0 );
      }
   }

   return(0);
}

//---------------------------------------------------------------------------------------------------------

String_t String_t::Get_Line (void)
{
   String_t s;

   if (str)
   {
      char *p  = str;
      char *pe = str + strlen(str);

      while ( (p<pe) && (*p!='\n') ) { p++; }

      int n = (p-str);

      if (*p=='\n') { p++; }

      if (n>0)
      {
         s.str = new char[n+1];

         char *p  =   str;
         char *q  = s.str;

         for (int i=0; (i<n); i++)  { q[i]= p[i]; }
         q[n]=0;

         p = str+n;
         q = str;

         if (*p=='\n') { p++; }
 
         while (p<pe) *q++ = *p++;

         *q++ = 0;
      }
   }

   return (s); 
}

//---------------------------------------------------------------------------------------------------------

char *String_t::Get_ENV (const char *s)
{
   if (str) { delete [] str; str=0; } 

   if (!s) { printf("String_t::Get_ENV() - null environment variable parameter"); return(0); }

   char *env = getenv(s);

   if (!env) { printf("String_t::Get_ENV() - environment varaible \"%s\" not set", s); return(0); }

   str = ( env ? strdup(env) : 0 );

   return(str);
}

//---------------------------------------------------------------------------------------------------------

int strlen(const String_t & s) { return( s.str ? strlen(s.str) : 0 ); }

