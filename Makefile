# Makefile for the ROOT test programs.
# This Makefile shows nicely how to compile and link applications
# using the ROOT libraries on all supported platforms.
#
# Copyright (c) 2000 Rene Brun and Fons Rademakers
#
# Author: Fons Rademakers, 29/2/2000

ARCH          = linux

CXX           =
ObjSuf        = o
SrcSuf        = cxx
ExeSuf        =
DllSuf        = so
OutPutOpt     = -o 

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)



ifeq ($(ARCH),linux)
# Linux with egcs, gcc 2.9x, gcc 3.x (>= RedHat 5.2)
CXX           = g++
CXXFLAGS      = -O -Wno-deprecated -fPIC -m32 -fno-inline -Wno-write-strings 
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -shared
INCLUDE	      = $(ROOTSYS)/include/

endif

CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS) $(SYSLIBS)
GLIBS         = $(ROOTGLIBS) $(SYSLIBS) 

#------------------------------------------------------------------------------
All : Fpdchan.so 


Fpdchan.so:	Fpdchan.h dataSet.h dataSet.cxx FpdMap.cxx CalibStr.cxx FpdchanLinkDef.h TObjFloat_t.h TObjFloat_t.cxx  Fill.cxx   FillData.cxx fOneD.h fOneD.cxx fEvt.h fEntry.h  CalibStr.h FpdMap.h Constants.h Fill.h FillData.h RunData.h RunData.cxx Filentry.h Filentry.cxx FilesSet.h FilesSet.cxx AnalTools.h AnalTools.cxx Geom.h Geom.cxx FitTower.cxx FitTower.h TowerUtil.cxx TowerUtil.h Trigger.h Trigger.cxx simh112.h simh112.cxx Sim.h Sim.cxx Cell.h Cell.cxx Gen.h Gen.cxx Qt.h Qt.cxx poutTree.h poutTree.cxx LVec.h LVec.cxx ClusterQt.h ClusterQt.cxx Yiqun.h Yiqun.cxx TrigQt.h TrigQt.cxx WasExternal.h WasExternal.cxx TowerFPD.h TowerFPD.cxx HitCluster.h HitCluster.cxx PhotonHitFPD.h PhotonHitFPD.cxx Vertex.h Vertex.cxx Sheet.h Sheet.cxx RunDepCor.h RunDepCor.cxx CellTDep.h CellTDep.cxx O2Output.h O2Output.cxx RunDepMgr.h RunDepMgr.cxx Asymmetry.h Asymmetry.cxx PullBBC.h PullBBC.cxx PullRP.h PullRP.cxx ProjectGCor.h ProjectGCor.cxx

		@echo "Generating dictionary $@..."
		@rootcint -f FpdchanDict.cxx -c Fpdchan.h CalibStr.h dataSet.h FpdMap.h  Fill.h FillData.h RunData.h Filentry.h FilesSet.h AnalTools.h Geom.h FitTower.h Trigger.h simh112.h Sim.h Cell.h Gen.h Qt.h poutTree.h LVec.h ClusterQt.h Yiqun.h TrigQt.h WasExternal.h TowerFPD.h PhotonHitFPD.h HitCluster.h TowerUtil.h Vertex.h Sheet.h  CellTDep.h RunDepCor.h O2Output.h RunDepMgr.h  Asymmetry.h PullBBC.h PullRP.h ProjectGCor.h FpdchanLinkDef.h 

		g++ dataSet.cxx FpdMap.cxx CalibStr.cxx FpdchanDict.cxx TObjFloat_t.cxx fOneD.cxx Fill.cxx FillData.cxx RunData.cxx  Filentry.cxx FilesSet.cxx AnalTools.cxx Geom.cxx FitTower.cxx TowerUtil.cxx Trigger.cxx simh112.cxx Sim.cxx Cell.cxx Gen.cxx Qt.cxx poutTree.cxx LVec.cxx ClusterQt.cxx Yiqun.cxx TrigQt.cxx WasExternal.cxx TowerFPD.cxx PhotonHitFPD.cxx HitCluster.cxx Vertex.cxx Sheet.cxx RunDepCor.cxx CellTDep.cxx O2Output.cxx RunDepMgr.cxx Asymmetry.cxx PullBBC.cxx PullRP.cxx ProjectGCor.cxx -shared -o Fpdchan.so $(CXXFLAGS) $(GLIBS) -lMinuit


