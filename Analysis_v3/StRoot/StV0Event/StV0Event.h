#ifndef StV0Event_h
#define StV0Event_h

#include "TMath.h"
#include "TObject.h"
#include "TClonesArray.h"
#include "StThreeVectorF.hh"
#include "TVector2.h"
#include "TLorentzVector.h"

class StV0Track : public TObject
{
  private:
    // Track properties
    // TrackA: for Lambda/anti-Lambda => proton/anti-proton
    // TrackB: for Lambda/anti-Lambda => pi_minus/pi_plus
    Float_t mMass2A; // mass2 of track A
    Float_t mMass2B;
    Float_t mNSigA; // nsigma dE/dx of particle A
    Float_t mNSigB;
    Float_t mDcaA; // distance of closest approach of particle A: dca*charge
    Float_t mDcaB;
    TLorentzVector mGTrackA; // Lorentz Vector for Global track A (px, py, pz, mass)
    TLorentzVector mGTrackB;
    TLorentzVector mPTrackA; // Lorentz Vector for Primary track A (px, py, pz, mass)
    TLorentzVector mPTrackB;
    Int_t mFlagA; // Flag for Event A: 0 for Same Event, others for Mixed Event
    Int_t mFlagB;
    Float_t mDcaAB; // dca between track A and track B
    Float_t mDecayLength; // decay length of mother particle
    Float_t mDcaV0; // dca between mother particle and primary vertex

  public:
    StV0Track() :
      mMass2A(-999.9),mMass2B(-999.9),mNSigA(-999.9),mNSigB(-999.9),mDcaA(-999.9),mDcaB(-999.9),mGTrackA(0,0,0,0),mGTrackB(0,0,0,0),mPTrackA(0,0,0,0),mPTrackB(0,0,0,0),mFlagA(-1),mFlagB(-1),mDcaAB(-999.9),mDecayLength(-999.9),mDcaV0(-999.9)
    {
    }
    ~StV0Track() {}

    // setters
    void setMass2A(Float_t f)                  { mMass2A = f;         }
    void setMass2B(Float_t f)                  { mMass2B = f;         }
    void setNSigA(Float_t f)                   { mNSigA = f;          }
    void setNSigB(Float_t f)                   { mNSigB = f;          }
    void setDcaA(Float_t f)                    { mDcaA = f;           }
    void setDcaB(Float_t f)                    { mDcaB = f;           }
    void setGTrackA(TLorentzVector f)          { mGTrackA = f;        }
    void setGTrackB(TLorentzVector f)          { mGTrackB = f;        }
    void setPTrackA(TLorentzVector f)          { mPTrackA = f;        }
    void setPTrackB(TLorentzVector f)          { mPTrackB = f;        }
    void setFlagA(Int_t f)                     { mFlagA = f;          }
    void setFlagB(Int_t f)                     { mFlagB = f;          }
    void setDcaAB(Float_t f)                   { mDcaAB = f;          }
    void setDecayLength(Float_t f)             { mDecayLength = f;    }
    void setDcaV0(Float_t f)                   { mDcaV0 = f;          }

    // getters
    Float_t getMass2A() const                  { return mMass2A;      }
    Float_t getMass2B() const                  { return mMass2B;      }
    Float_t getNSigA() const                   { return mNSigA;       }
    Float_t getNSigB() const                   { return mNSigB;       }
    Float_t getDcaA() const                    { return mDcaA;        }
    Float_t getDcaB() const                    { return mDcaB;        }
    TLorentzVector getGTrackA() const          { return mGTrackA;     }
    TLorentzVector getGTrackB() const          { return mGTrackB;     }
    TLorentzVector getPTrackA() const          { return mPTrackA;     }
    TLorentzVector getPTrackB() const          { return mPTrackB;     }
    Int_t getFlagA() const                     { return mFlagA;       }
    Int_t getFlagB() const                     { return mFlagB;       }
    Float_t getDcaAB() const                   { return mDcaAB;       }
    Float_t getDecayLength() const             { return mDecayLength; }
    Float_t getDcaV0() const                   { return mDcaV0;       }

    ClassDef(StV0Track,1)  // A simple track of a particle
};

class StV0Event: public TObject
{
  private:
    StThreeVectorF mPrimaryvertex;
    Int_t   mRunId;
    Int_t   mRefMult;
    Int_t   mCentrality;
    Int_t   mN_prim;
    Int_t   mN_non_prim;
    Int_t   mN_Tof_match;

    Float_t mZDCx;
    Float_t mBBCx;
    Float_t mVzVpd;

    UShort_t      fNumTracks;

    TVector2 mQ2East[4];
    TVector2 mQ2West[4];
    TVector2 mQ3East[4];
    TVector2 mQ3West[4];
    Int_t   mNumTrackEast[4];
    Int_t   mNumTrackWest[4];

    TClonesArray* fTracks;      //->

  public:
    StV0Event() :
      mPrimaryvertex(-1.0,-1.0,-1.0),mRunId(-1),mRefMult(-1),mCentrality(-1),mN_prim(-1),mN_non_prim(-1),mN_Tof_match(-1),mZDCx(-1),mBBCx(-1),mVzVpd(-1),fNumTracks(0)
    {
      for(Int_t j = 0; j < 4; j++)
      {
	mQ2East[j].Set(-999.9,-999.9); // QVector2 East
	mQ2West[j].Set(-999.9,-999.9); // QVector2 West
	mQ3East[j].Set(-999.9,-999.9); // QVector3 East
	mQ3West[j].Set(-999.9,-999.9); // QVector3 West

	mNumTrackEast[j] = 0;
	mNumTrackWest[j] = 0;
      }
      fTracks      = new TClonesArray( "StV0Track", 10 );
    }
    ~StV0Event()
    {
      delete fTracks;
      fTracks = NULL;
    }

    void       setPrimaryVertex(StThreeVectorF r)      { mPrimaryvertex = r;     }
    StThreeVectorF    getPrimaryVertex() const         { return mPrimaryvertex;  }

    void       setRunId(Int_t  r)                      { mRunId = r;             }
    Int_t      getRunId() const                        { return mRunId;          }

    void       setRefMult(Int_t r)                     { mRefMult = r;           }
    Int_t      getRefMult() const                      { return mRefMult;        }

    void       setCentrality(Int_t r)                  { mCentrality = r;        }
    Int_t      getCentrality() const                   { return mCentrality;     }

    void       setN_prim(Int_t r)                      { mN_prim = r;            }
    Int_t      getN_prim() const                       { return mN_prim;         }

    void       setN_non_prim(Int_t r)                  { mN_non_prim = r;        }
    Int_t      getN_non_prim() const                   { return mN_non_prim;     }

    void       setN_Tof_match(Int_t r)                 { mN_Tof_match = r;       }
    Int_t      getN_Tof_match() const                  { return mN_Tof_match;    }


    void       setZDCx(Float_t r)                      { mZDCx = r;              }
    Float_t    getZDCx() const                         { return mZDCx;           }

    void       setBBCx(Float_t r)                      { mBBCx = r;              }
    Float_t    getBBCx() const                         { return mBBCx;           }

    void       setVzVpd(Float_t r)                     { mVzVpd = r;             }
    Float_t    getVzVpd() const                        { return mVzVpd;          }

    // ---------------------------------------QVector---------------------------------------------
    // QVector2 East
    void       setQ2East(TVector2 r, Int_t j)          { mQ2East[j] = r;         }
    TVector2   getQ2East(Int_t j) const                { return mQ2East[j];      }
    // QVector2 West
    void       setQ2West(TVector2 r, Int_t j)          { mQ2West[j] = r;         }
    TVector2   getQ2West(Int_t j) const                { return mQ2West[j];      }
    // QVector3 East
    void       setQ3East(TVector2 r, Int_t j)          { mQ3East[j] = r;         }
    TVector2   getQ3East(Int_t j) const                { return mQ3East[j];      }
    // QVector3 West
    void       setQ3West(TVector2 r, Int_t j)          { mQ3West[j] = r;         }
    TVector2   getQ3West(Int_t j) const                { return mQ3West[j];      }
    // ---------------------------------------QVector---------------------------------------------

    // -----------------------------------Number of Tracks----------------------------------------
    // East
    void       setNumTrackEast(Int_t r, Int_t j)       { mNumTrackEast[j] = r;   }
    Int_t      getNumTrackEast(Int_t j) const          { return mNumTrackEast[j];}
    // West
    void       setNumTrackWest(Int_t r, Int_t j)       { mNumTrackWest[j] = r;   }
    Int_t      getNumTrackWest(Int_t j) const          { return mNumTrackWest[j];}
    // -----------------------------------Number of Tracks----------------------------------------
    StV0Track* createTrack()
    {
      if (fNumTracks == fTracks->GetSize())
	fTracks->Expand( fNumTracks + 10 );
      if (fNumTracks >= 10000)
      {
	Fatal( "StV0Event::createTrack()", "ERROR: Too many tracks (>10000)!" );
	exit( 2 );
      }

      new((*fTracks)[fNumTracks++]) StV0Track;
      return (StV0Track*)((*fTracks)[fNumTracks - 1]);
    }
    void clearTrackList()
    {
      fNumTracks   = 0;
      fTracks      ->Clear();
    }
    UShort_t getNumTracks() const
    {
      return fNumTracks;
    }
    StV0Track* getTrack(UShort_t i) const
    {
      return i < fNumTracks ? (StV0Track*)((*fTracks)[i]) : NULL;
    }

    ClassDef(StV0Event,1)  // A simple event compiled of tracks
};


#endif
