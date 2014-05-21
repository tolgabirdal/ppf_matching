
#pragma once

#include "opencv2/core.hpp"

class RGBDFeatureDetector : public virtual Algorithm
{
public:
	RGBDMatcher(void);
	~RGBDMatcher(void);
};

class CV_EXPORTS_W_SIMPLE FeaturePoint3D
{
public:
    //! the default constructor
    CV_WRAP FeaturePoint3D() : pt(0,0, 0), size(0), angle(-1), response(0), octave(0), class_id(-1) {}
    //! the full constructor
    KeyPoint(Point2f _pt, float _size, float _angle=-1,
            float _response=0, int _octave=0, int _class_id=-1)
            : pt(_pt), size(_size), angle(_angle),
            response(_response), octave(_octave), class_id(_class_id) {}
    //! another form of the full constructor
    CV_WRAP KeyPoint(float x, float y, float _size, float _angle=-1,
            float _response=0, int _octave=0, int _class_id=-1)
            : pt(x, y), size(_size), angle(_angle),
            response(_response), octave(_octave), class_id(_class_id) {}

    size_t hash() const;

    CV_PROP_RW Point2f pt; //!< coordinates of the keypoints
    CV_PROP_RW float size; //!< diameter of the meaningful keypoint neighborhood
    CV_PROP_RW float angle; //!< computed orientation of the keypoint (-1 if not applicable);
                            //!< it's in [0,360) degrees and measured relative to
                            //!< image coordinate system, ie in clockwise.
    CV_PROP_RW float response; //!< the response by which the most strong keypoints have been selected. Can be used for the further sorting or subsampling
    CV_PROP_RW int octave; //!< octave (pyramid layer) from which the keypoint has been extracted
    CV_PROP_RW int class_id; //!< object class (if the keypoints need to be clustered by an object they belong to)
};

class CV_EXPORTS_W PPFFeatureDetector : public RGBDFeatureDetector
{
public:
    /*
     * Detect keypoints in an image.
     * image        The image.
     * keypoints    The detected keypoints.
     * mask         Mask specifying where to look for keypoints (optional). Must be a char
     *              matrix with non-zero values in the region of interest.
     * useProvidedKeypoints If true, the method will skip the detection phase and will compute
     *                      descriptors for the provided keypoints
     */
    CV_WRAP_AS(detectAndCompute) virtual void operator()( InputArray image, CV_OUT vector<KeyPoint>& keypoints,
                                     OutputArray descriptors,
                                     bool useProvidedKeypoints=false ) const = 0;

    // Create feature detector and descriptor extractor by name.
    CV_WRAP static Ptr<Feature2D> create( const string& name );
};