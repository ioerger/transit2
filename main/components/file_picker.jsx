/** @jsx html */
/// <reference no-default-lib="true"/>
/// <reference lib="dom" />
/// <reference lib="dom.asynciterable" />
/// <reference lib="deno.ns" />
import { html } from "../elemental.js"
import { Column, Row, Input, Code, EasyFilePicker, askForFiles } from "../elements.jsx"


export const EasyFilePicker = ({children, onChange, defaultWidth="16rem", backgroundColor="#939393", ...otherArgs})=>{
    let codeElement
    return <Row 
        style={`
            --default-width: ${defaultWidth};
            width: var(--default-width);
            max-width: var(--default-width);
            overflow-x: visible;
        `}
        >
            {codeElement = <Code
                style={`
                    width: var(--default-width);
                    max-width: var(--default-width);
                    display: block;
                    overflow-x: auto;
                    background: ${backgroundColor};
                    position: relative;
                    box-shadow: 0 4px 5px 0 rgba(0,0,0,0), 0 1px 10px 0 rgba(0,0,0,0), 0 2px 4px -1px rgba(0,0,0,0);
                    margin-bottom: 0.5rem;
                `}
                hoverStyle={`
                    min-width: max-content;
                    overflow-x: visible;
                    box-shadow: 0 4px 5px 0 rgba(0,0,0,0.14), 0 1px 10px 0 rgba(0,0,0,0.12), 0 2px 4px -1px rgba(0,0,0,0.3);
                `}
                onClick={async (event)=>{
                    const files = await askForFiles()
                    if (files && files.length > 0) {
                        codeElement.innerText = files[0].name
                        onChange(files[0])
                    }
                }}
                {...otherArgs}
                >
                    <Row horizontalAlignment="center" width="100%" >
                        {children}
                    </Row>
            </Code>}
    </Row>
}