module.exports = {
  "id": "venn.js",
  "viewports": [
    {
      "label": "phone_portrait",
      "width": 320,
      "height": 480
    },
    {
      "label": "phone_landscaope",
      "width": 480,
      "height": 320
    },
    {
      "label": "tablet",
      "width": 1024,
      "height": 768
    }
  ],
  "onBeforeScript": "puppet/onBefore.js",
  "onReadyScript": "puppet/onReady.js",
  "scenarios": [
    {
      "label": "simple",
      "url": "./examples/simple.html",
      "delay": 500,
      "selectors": ["viewport"],
      "misMatchThreshold" : 0.1,
      "requireSameDimensions": true
    },
    {
      "label": "simple_viewbox",
      "url": "./examples/simple_viewbox.html",
      "delay": 500,
      "selectors": ["viewport"],
      "misMatchThreshold" : 0.1,
      "requireSameDimensions": true
    }
  ],
  "paths": {
    "bitmaps_reference": "backstop_data/bitmaps_reference",
    "bitmaps_test": "backstop_data/bitmaps_test",
    "engine_scripts": "backstop_data/engine_scripts",
    "html_report": "backstop_data/html_report",
    "ci_report": "backstop_data/ci_report"
  },
  "report": ["browser"],
  "engine": "puppeteer",
  "engineOptions": {
    "args": ["--no-sandbox"]
  },
  "asyncCaptureLimit": 5,
  "asyncCompareLimit": 50,
  "debug": false,
  "debugWindow": false
};
