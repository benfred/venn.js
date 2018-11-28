import resolve from 'rollup-plugin-node-resolve';
import commonjs from 'rollup-plugin-commonjs';

export default {
    input: 'index.js',
    output: {
        file: 'build/venn.js',
        format: 'umd',
        name: 'venn',
        globals: {
            'd3-selection': 'd3',
            'd3-transition': 'd3'
        }
    },
    plugins: [
        resolve({
            jsnext: true,
            only: ['fmin'],
            main: true
        }),
        commonjs()
    ]
};
