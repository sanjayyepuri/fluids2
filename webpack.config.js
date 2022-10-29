const path = require("path");
const HtmlWebpackPlugin = require("html-webpack-plugin");

module.exports = {
    entry: "./src/index.js",
    output: {
        path: path.resolve(__dirname, "dist"),
        clean: true,
    },
    plugins:[
        new HtmlWebpackPlugin({title: "Fluids2", template: "src/index.html"})
    ],
};