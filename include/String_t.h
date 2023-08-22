#pragma once

#ifndef __PDS_STRING_H
#define __PDS_STRING_H

//---------------------------------------------------------------------------------------------------------

class String_t
{
   public:
      char *str;

   public:
       String_t(void);
       String_t(const char *s0  , const char *s1=0, const char *s2=0, const char *s3=0, const char *s4=0,
                const char *s5=0, const char *s6=0, const char *s7=0, const char *s8=0, const char *s9=0 );
       String_t(const String_t & s);

      ~String_t();

      const String_t & operator = (const String_t & s);
      const String_t & operator = (const char      *s);

      int   Length() const;

      operator          char     *() const;
      operator          int       () const;
      operator unsigned int       () const;
      operator          short     () const;
      operator unsigned short     () const;
      operator          long      () const;
      operator unsigned long      () const;
      operator          long long () const;
      operator unsigned long long () const;
      operator          float     () const;
      operator          double    () const;

      void       operator += (const char    *s);
      void       operator += (const String_t  &s);

      int        operator == (const char      *s) const;
      int        operator != (const char      *s) const;
      int        operator <  (const char      *s) const;
      int        operator >  (const char      *s) const;
      int        operator == (const String_t  &s) const;
      int        operator != (const String_t  &s) const;
      int        operator <  (const String_t  &s) const;
      int        operator >  (const String_t  &s) const;

      int        Equiv       (const char *s) const;

      void       Append      (const char *s0  , const char *s1=0, const char *s2=0, const char *s3=0, const char *s4=0,
                              const char *s5=0, const char *s6=0, const char *s7=0, const char *s8=0, const char *s9=0 );

      void       Upper       (void);
      void       Lower       (void);

      void       Replace     (const char c0, const char c1);

      void       Strip       (void);
      void       Strip_L     (void);
      void       Strip_R     (void);

      int        Prefix      (const char *pfx) const;
      int        Suffix      (const char *sfx) const;

      void       Extract     (const int i0,               String_t & s);
      void       Extract     (const int i0, const int i1, String_t & s);

      char      *Find        (const char  c) const;
      char      *Find        (const char *s) const;

      String_t   Get_Line    (void);

      char      *Get_ENV     (const char *s);

};

extern int strlen(const String_t & s);

#endif  // __PDS_STRING_H
