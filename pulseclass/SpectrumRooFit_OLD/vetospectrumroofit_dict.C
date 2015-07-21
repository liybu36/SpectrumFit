// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME vetospectrumroofit_dict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "../vetopulsesplit.hh"
#include "../odselector/odselector.h"
#include "../DSTtreeSelector/DSTtreeSelector.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_vetopulsesplit(void *p = 0);
   static void *newArray_vetopulsesplit(Long_t size, void *p);
   static void delete_vetopulsesplit(void *p);
   static void deleteArray_vetopulsesplit(void *p);
   static void destruct_vetopulsesplit(void *p);
   static void streamer_vetopulsesplit(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::vetopulsesplit*)
   {
      ::vetopulsesplit *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::vetopulsesplit >(0);
      static ::ROOT::TGenericClassInfo 
         instance("vetopulsesplit", ::vetopulsesplit::Class_Version(), "../vetopulsesplit.hh", 16,
                  typeid(::vetopulsesplit), DefineBehavior(ptr, ptr),
                  &::vetopulsesplit::Dictionary, isa_proxy, 16,
                  sizeof(::vetopulsesplit) );
      instance.SetNew(&new_vetopulsesplit);
      instance.SetNewArray(&newArray_vetopulsesplit);
      instance.SetDelete(&delete_vetopulsesplit);
      instance.SetDeleteArray(&deleteArray_vetopulsesplit);
      instance.SetDestructor(&destruct_vetopulsesplit);
      instance.SetStreamerFunc(&streamer_vetopulsesplit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::vetopulsesplit*)
   {
      return GenerateInitInstanceLocal((::vetopulsesplit*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::vetopulsesplit*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_odselector(void *p = 0);
   static void *newArray_odselector(Long_t size, void *p);
   static void delete_odselector(void *p);
   static void deleteArray_odselector(void *p);
   static void destruct_odselector(void *p);
   static void streamer_odselector(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::odselector*)
   {
      ::odselector *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::odselector >(0);
      static ::ROOT::TGenericClassInfo 
         instance("odselector", ::odselector::Class_Version(), "../odselector/odselector.h", 25,
                  typeid(::odselector), DefineBehavior(ptr, ptr),
                  &::odselector::Dictionary, isa_proxy, 16,
                  sizeof(::odselector) );
      instance.SetNew(&new_odselector);
      instance.SetNewArray(&newArray_odselector);
      instance.SetDelete(&delete_odselector);
      instance.SetDeleteArray(&deleteArray_odselector);
      instance.SetDestructor(&destruct_odselector);
      instance.SetStreamerFunc(&streamer_odselector);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::odselector*)
   {
      return GenerateInitInstanceLocal((::odselector*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::odselector*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_DSTtreeSelector(void *p = 0);
   static void *newArray_DSTtreeSelector(Long_t size, void *p);
   static void delete_DSTtreeSelector(void *p);
   static void deleteArray_DSTtreeSelector(void *p);
   static void destruct_DSTtreeSelector(void *p);
   static void streamer_DSTtreeSelector(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::DSTtreeSelector*)
   {
      ::DSTtreeSelector *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::DSTtreeSelector >(0);
      static ::ROOT::TGenericClassInfo 
         instance("DSTtreeSelector", ::DSTtreeSelector::Class_Version(), "../DSTtreeSelector/DSTtreeSelector.h", 25,
                  typeid(::DSTtreeSelector), DefineBehavior(ptr, ptr),
                  &::DSTtreeSelector::Dictionary, isa_proxy, 16,
                  sizeof(::DSTtreeSelector) );
      instance.SetNew(&new_DSTtreeSelector);
      instance.SetNewArray(&newArray_DSTtreeSelector);
      instance.SetDelete(&delete_DSTtreeSelector);
      instance.SetDeleteArray(&deleteArray_DSTtreeSelector);
      instance.SetDestructor(&destruct_DSTtreeSelector);
      instance.SetStreamerFunc(&streamer_DSTtreeSelector);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::DSTtreeSelector*)
   {
      return GenerateInitInstanceLocal((::DSTtreeSelector*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::DSTtreeSelector*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr vetopulsesplit::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *vetopulsesplit::Class_Name()
{
   return "vetopulsesplit";
}

//______________________________________________________________________________
const char *vetopulsesplit::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::vetopulsesplit*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int vetopulsesplit::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::vetopulsesplit*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *vetopulsesplit::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::vetopulsesplit*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *vetopulsesplit::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::vetopulsesplit*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr odselector::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *odselector::Class_Name()
{
   return "odselector";
}

//______________________________________________________________________________
const char *odselector::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::odselector*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int odselector::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::odselector*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *odselector::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::odselector*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *odselector::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::odselector*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr DSTtreeSelector::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *DSTtreeSelector::Class_Name()
{
   return "DSTtreeSelector";
}

//______________________________________________________________________________
const char *DSTtreeSelector::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::DSTtreeSelector*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int DSTtreeSelector::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::DSTtreeSelector*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *DSTtreeSelector::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::DSTtreeSelector*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *DSTtreeSelector::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::DSTtreeSelector*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void vetopulsesplit::Streamer(TBuffer &R__b)
{
   // Stream an object of class vetopulsesplit.

   TObject::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_vetopulsesplit(void *p) {
      return  p ? new(p) ::vetopulsesplit : new ::vetopulsesplit;
   }
   static void *newArray_vetopulsesplit(Long_t nElements, void *p) {
      return p ? new(p) ::vetopulsesplit[nElements] : new ::vetopulsesplit[nElements];
   }
   // Wrapper around operator delete
   static void delete_vetopulsesplit(void *p) {
      delete ((::vetopulsesplit*)p);
   }
   static void deleteArray_vetopulsesplit(void *p) {
      delete [] ((::vetopulsesplit*)p);
   }
   static void destruct_vetopulsesplit(void *p) {
      typedef ::vetopulsesplit current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_vetopulsesplit(TBuffer &buf, void *obj) {
      ((::vetopulsesplit*)obj)->::vetopulsesplit::Streamer(buf);
   }
} // end of namespace ROOT for class ::vetopulsesplit

//______________________________________________________________________________
void odselector::Streamer(TBuffer &R__b)
{
   // Stream an object of class odselector.

   TSelector::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_odselector(void *p) {
      return  p ? new(p) ::odselector : new ::odselector;
   }
   static void *newArray_odselector(Long_t nElements, void *p) {
      return p ? new(p) ::odselector[nElements] : new ::odselector[nElements];
   }
   // Wrapper around operator delete
   static void delete_odselector(void *p) {
      delete ((::odselector*)p);
   }
   static void deleteArray_odselector(void *p) {
      delete [] ((::odselector*)p);
   }
   static void destruct_odselector(void *p) {
      typedef ::odselector current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_odselector(TBuffer &buf, void *obj) {
      ((::odselector*)obj)->::odselector::Streamer(buf);
   }
} // end of namespace ROOT for class ::odselector

//______________________________________________________________________________
void DSTtreeSelector::Streamer(TBuffer &R__b)
{
   // Stream an object of class DSTtreeSelector.

   TSelector::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_DSTtreeSelector(void *p) {
      return  p ? new(p) ::DSTtreeSelector : new ::DSTtreeSelector;
   }
   static void *newArray_DSTtreeSelector(Long_t nElements, void *p) {
      return p ? new(p) ::DSTtreeSelector[nElements] : new ::DSTtreeSelector[nElements];
   }
   // Wrapper around operator delete
   static void delete_DSTtreeSelector(void *p) {
      delete ((::DSTtreeSelector*)p);
   }
   static void deleteArray_DSTtreeSelector(void *p) {
      delete [] ((::DSTtreeSelector*)p);
   }
   static void destruct_DSTtreeSelector(void *p) {
      typedef ::DSTtreeSelector current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_DSTtreeSelector(TBuffer &buf, void *obj) {
      ((::DSTtreeSelector*)obj)->::DSTtreeSelector::Streamer(buf);
   }
} // end of namespace ROOT for class ::DSTtreeSelector

namespace {
  void TriggerDictionaryInitialization_vetospectrumroofit_dict_Impl() {
    static const char* headers[] = {
"../vetopulsesplit.hh",
"../odselector/odselector.h",
"../DSTtreeSelector/DSTtreeSelector.h",
0
    };
    static const char* includePaths[] = {
"/Users/hqian/Documents/root6/root-6.02.05/include",
"/Users/hqian/OneDrive/montecarlo/g4ds10/Linux-g++/pulse_splitter/pulseclass/SpectrumRooFit/",
0
    };
    static const char* fwdDeclCode = 
R"DICTFWDDCLS(
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$../vetopulsesplit.hh")))  vetopulsesplit;
class __attribute__((annotate("$clingAutoload$../odselector/odselector.h")))  odselector;
class __attribute__((annotate("$clingAutoload$../DSTtreeSelector/DSTtreeSelector.h")))  DSTtreeSelector;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "../vetopulsesplit.hh"
#include "../odselector/odselector.h"
#include "../DSTtreeSelector/DSTtreeSelector.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"DSTtreeSelector", payloadCode, "@",
"odselector", payloadCode, "@",
"vetopulsesplit", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("vetospectrumroofit_dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_vetospectrumroofit_dict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_vetospectrumroofit_dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_vetospectrumroofit_dict() {
  TriggerDictionaryInitialization_vetospectrumroofit_dict_Impl();
}
