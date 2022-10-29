const path = require("path");
const HtmlWebpackPlugin = require("html-webpack-plugin");

module.exports = {
    entry: "./src/index.ts",
    module: {
        rules: [
            {
                test: /\.tsx?$/,
                use: "ts-loader",
                exclude: /node_modules/, 
            },
        ],
    },
    resolve: {
        extensions: ['.tsx', '.ts', '.js'],
    },
    output: {
        path: path.resolve(__dirname, "dist"),
        clean: true,
    },
    plugins:[
        new HtmlWebpackPlugin({title: "Fluids2", template: "./index.html"})
    ],
    experiments: {
        asyncWebAssembly: true,
    },
    mode: "development",
};