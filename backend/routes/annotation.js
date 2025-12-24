const express = require('express');
const router = express.Router();
const { handleAnnotation, handleAnnotationFromPrediction } = require('../services/annotationService');

router.get('/api/run-annotation', (req, res) => handleAnnotation(req, res, 'GET'));
router.post('/api/run-annotation', (req, res) => handleAnnotation(req, res, 'POST'));

router.get('/api/run-annotation-from-prediction', (req, res) =>
  handleAnnotationFromPrediction(req, res, 'GET')
);
router.post('/api/run-annotation-from-prediction', (req, res) =>
  handleAnnotationFromPrediction(req, res, 'POST')
);

module.exports = router;