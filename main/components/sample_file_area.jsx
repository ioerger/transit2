/** @jsx html */
/// <reference no-default-lib="true"/>
/// <reference lib="dom" />
/// <reference lib="dom.asynciterable" />
/// <reference lib="deno.ns" />
import { html } from "../elemental.js"
import { Column, Row, Input, askForFiles } from "../elements.jsx"


export const SampleFilePicker = ({ children, style, }) => {
    return <Column padding="3rem">
        <button
            onClick={async (event)=>{
                const comwigFiles = askForFiles()
                for (const each of comwigFiles) {
                    askForFiles()
                }
                console.log("file selected", )
            }} 
            margin="2rem"
            >
                Select File
        </button>
        <Column id="wigTable" background="gray" padding="3rem">
            Table
        </Column>
    </Column>
}

