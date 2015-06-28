module.exports = function(grunt) {
    grunt.initConfig({
        pkg: grunt.file.readJSON('package.json'),

        concat: {
            dist: {
                src: [
                    'src/init.js',
                    'src/diagram.js',
                    'src/layout.js',
                    'src/fmin.js',
                    'src/circleintersection.js',
                    'src/export.js',
                ],
                dest: 'venn.js',
            }
        },

        uglify: {
            build: {
                src: 'venn.js',
                dest: 'venn.min.js'
            }
        },

        qunit: {
            all: ['tests/unittest.html'],
        },

        jshint: {
            all: ['Gruntfile.js', 'src/*.js', 'tests/*js'],
        }
    });

    grunt.loadNpmTasks('grunt-contrib-concat');
    grunt.loadNpmTasks('grunt-contrib-uglify');
    grunt.loadNpmTasks('grunt-contrib-qunit');
    grunt.loadNpmTasks('grunt-contrib-jshint');
    grunt.registerTask('default', ['concat', 'uglify', 'jshint', 'qunit']);
};
