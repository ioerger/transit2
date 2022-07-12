/** @jsx html */
/// <reference no-default-lib="true"/>
/// <reference lib="dom" />
/// <reference lib="dom.asynciterable" />
/// <reference lib="deno.ns" />
import { html } from "../elemental.js"
import { Event, trigger, everyTime, once } from "https://deno.land/x/good@0.5.15/events.js"
import { Column, Row, Input, Code, askForFiles } from "../elements.jsx"
import { EasyFilePicker } from "./file_picker.jsx"
import { data, events } from "../data.jsx"

export const Annotation = ({ children, style, }) => {
    return <Column>
        <span class="custom-header">
            Annotation File
        </span>
        <Row
            padding="1.2rem 1rem"
            horizontalAlignment={`space-between`}
            verticalAlignment={"center"}
            background="whitesmoke"
            border="lightgray solid 1px"
            min-width="fit-content"
            border-radius="0.7rem"
            >
                <EasyFilePicker
                    defaultWidth="18rem"
                    onChange={(file)=>{
                        data.annotation = file.name
                        console.log(`data.annotation is:`,data.annotation)
                        trigger(events.annotationAdded)
                    }}
                    >
                        [Click to add Annotation File]
                </EasyFilePicker>
        </Row>
    </Column>
}

