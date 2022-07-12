/** @jsx html */
/// <reference no-default-lib="true"/>
/// <reference lib="dom" />
/// <reference lib="dom.asynciterable" />
/// <reference lib="deno.ns" />
import { html, css, Elemental } from "https://raw.githubusercontent.com/jeff-hykin/elemental/f9764dd95cb0e645e1dffbe369b7a6467d2e77e3/main/deno.js"
import { capitalize, indent, toCamelCase, numberToEnglishArray, toPascalCase, toKebabCase, toSnakeCase, toScreamingtoKebabCase, toScreamingtoSnakeCase, toRepresentation, toString } from "https://deno.land/x/good@0.5.15/string.js"

// roadmap of tools:
    // Code
    // Collaspable
    // NestedMenus
    // Table
    // inputs:
        // Button
        // Dropdown
        // Search
            // MultiSelect
            // SingleSelect
        // ExpandingTextbox
        // Checkbox
        // PickOne
        // PickMany
        // Number
        // Slider
        // DatePicker
        // TimePicker
        // DateTimePicker
        // PhoneNumber
        // Email
        // Address
    // outputs:
        // Toast
        // Popover
        // Chip
        // LoadingSpinner
            // progress option
        // Tabs
    // Video
        // mp4/avi source link
    // ContextMenu
    

window.Elemental = Elemental // for debugging only

const randomId = Elemental.randomId
const combineClasses = Elemental.combineClasses

const translateAlignment = (name) => {
    if (name == "top" || name == "left") {
        return "flex-start"
    } else if (name == "bottom" || name == "right") {
        return "flex-end"
    } else {
        return name
    }
}

const classIds = {
    column: randomId(`column`),
    row: randomId(`row`),
    popUp: randomId(`popUp`),
    button: randomId(`button`),
    code: randomId(`code`),
}

document.body.appendChild(<style>{`
    .${classIds.column} {
        display: flex;
        flex-direction: column;
        transition: all 0.5s ease-in-out 0s;
    }
    .${classIds.row} {
        display: flex;
        flex-direction: row;
        transition: all 0.5s ease-in-out 0s;
    }
    .${classIds.popUp} {
        align-items: center;
        justify-content: center;
        position: absolute;
        top:0;
        left: 0;
        z-index: 5000;
        display: flex;
        width: 100vw;
        height: 100vh;
    }
    .${classIds.button} {
        background: whitesmoke;
        position: relative;
        box-shadow: 0 4px 5px 0 rgba(0,0,0,0.14), 0 1px 10px 0 rgba(0,0,0,0.12), 0 2px 4px -1px rgba(0,0,0,0.3);
        transition: all 0.5s ease-in-out 0s;
    }
    .${classIds.code} {
        white-space: pre;
        background: gray;
        position: relative;
        color: white;
        font-size: 11pt;
        padding: 0.5rem 0.7rem;
        border-radius: 0.4rem;
        transition: all 0.5s ease-in-out 0s;
    }
`}</style>)


const hoverStyleHelper = ({ element, hoverStyle }) => {
    if (hoverStyle) {
        let hoverStyleAlreadyActive = false
        const helper = document.createElement("div")
        const hoverStyleAsString = `${css(hoverStyle)}`
        helper.style.cssText = hoverStyleAsString // style string values change when attached to actual elements
        const styleObject = {}
        const keys = Object.values(helper.style) // yes I know it says keys= .values() but its true
        for (const key of keys) {
            styleObject[key] = helper.style[key]
        }
        const valuesBefore = {}
        
        element.addEventListener("mouseover", ()=>{
            if (!hoverStyleAlreadyActive) {
                hoverStyleAlreadyActive = true
                for (const key of keys) {
                    valuesBefore[key] = element.style[key]
                }
                element.style.cssText += hoverStyleAsString
            }
        })

        element.addEventListener("mouseout", ()=>{
            if (hoverStyleAlreadyActive) {
                hoverStyleAlreadyActive = false
                const style = element.style
                const mixinStyleObject = {}
                for (const [key, value] of Object.entries(styleObject)) {
                    // if it wasn't changed
                    if (style[key] == value) {
                        // then restore the old value
                        mixinStyleObject[key] = valuesBefore[key]
                    }
                }
                const mixinStyles = `${css(mixinStyleObject)}`
                style.cssText += mixinStyles           // this is needed for values with !important
                Object.assign(style, mixinStyleObject) // needed for empty values 
            }
        })
    }
}


export function Box({
        children,
        style,
        hoverStyle,
        row,
        column=true,
        center=false,
        verticalAlignment=null,
        horizontalAlignment=null,
        onBlur,
        onChange,
        onClick,
        onContextMenu,
        onDblClick,
        onMouseDown,
        onMouseEnter,
        onMouseLeave,
        onMouseMove,
        onMouseOut,
        onMouseOver,
        onMouseUp,
        ...otherArgs
    }) {
        let justify, align, text, theClass
        
        if (!row) {
            if (center) {
                verticalAlignment = verticalAlignment || "center"
                horizontalAlignment = horizontalAlignment || "center"
            } else {
                verticalAlignment = verticalAlignment || "top"
                horizontalAlignment = horizontalAlignment || "left"
            }
            justify = verticalAlignment
            align = horizontalAlignment
            text = horizontalAlignment
            theClass = classIds.column
        } else {
            if (center) {
                verticalAlignment = verticalAlignment || "center"
                horizontalAlignment = horizontalAlignment || "center"
            } else {
                verticalAlignment = verticalAlignment || "top"
                horizontalAlignment = horizontalAlignment || "left"
            }
            justify = horizontalAlignment
            align = verticalAlignment
            text = horizontalAlignment
            theClass = classIds.row
        }
        
        const element = <div
            class={combineClasses(theClass, otherArgs.class)}
            style={`justify-content: ${translateAlignment(justify)}; align-items: ${translateAlignment(align)}; text-align: ${text}; ${css(style)}; ${css(otherArgs)};`}
            onBlur={onBlur}
            onChange={onChange}
            onClick={onClick}
            onContextMenu={onContextMenu}
            onDblClick={onDblClick}
            onMouseDown={onMouseDown}
            onMouseEnter={onMouseEnter}
            onMouseLeave={onMouseLeave}
            onMouseMove={onMouseMove}
            onMouseOut={onMouseOut}
            onMouseOver={onMouseOver}
            onMouseUp={onMouseUp}
            {...otherArgs}
            >
                {children}
        </div>

        hoverStyleHelper({ element, hoverStyle })
        
        return element
}

export const Column = Box

export const Row = (arg)=>Box({...arg, row: true})

export const Input = ({ children, style, checked, value, ...otherArgs }) => {
    const element = <input
        class={otherArgs.class}
        style={`${css(style)}; ${css(otherArgs)};`}
        {...otherArgs}
        />
    // are not just html attributes
    Object.assign(element, {
        checked,
        value,
    })
    return element
}

export const Code = ({ children, style, hoverStyle, onMouseOver, onMouseOut, onClick, ...otherArgs }) => {
    const element = <code
        class={combineClasses(classIds.code, otherArgs.class)}
        style={`${css(style)}; ${css(otherArgs)};`}
        onClick={onClick}
        onMouseOver={onMouseOver}
        onMouseOut={onMouseOut}
        {...otherArgs}
        >
            {children}
    </code>

    hoverStyleHelper({ element, hoverStyle })

    return element
}

export const EasyFilePicker = ({ children, style, hoverStyle, onChange, ...otherArgs }) => {
    let element
    return element = <Code
        style={style}
        hoverStyle={hoverStyle}
        onClick={async (event)=>{
            const files = await askForFiles()
            if (files && files.length > 0) {
                element.innerText = files[0].name
                onChange(files[0])
            }
        }}
        {...otherArgs}
        >
            {children}
    </Code>
}

export const popUp = async ({ children, style, center, ...otherArgs })=>{
    const container = <div
        class={combineClasses(classIds.popUp, otherArgs.class)}
        onClick={event=>{
            // if actually clicked the container itself
            if (event.target == container) {
                // close the popUp
                container.remove()
            }
        }}
        >
            <Column verticalAlignment="top" horizontalAlignment="center" style="width: fit-content; height: 50vh; overflow-y: auto;">
                {children}
            </Column>
    </div>
    document.body.prepend(container)
    return container
}

export const askForFiles = async ()=>{
    return new Promise((resolve, reject)=>{
        const cleanResolve = (returnValue)=>{
            resolve(returnValue)
            window.removeEventListener("focus", listener)
            document.body.removeChild(filePicker)
        }
        const listener = ()=>cleanResolve([])
        window.addEventListener("focus", listener)
        let filePicker = <input
            type="file"
            onInput={event=>{ cleanResolve(event.target.files) }}
            onBlur={event=>{ cleanResolve([]) }}
            hidden
            />
        document.body.appendChild(filePicker)
        filePicker.click()
    })
}

export function Button({
        children,
        style,
        hoverStyle,
        row,
        column=true,
        center=false,
        verticalAlignment=null,
        horizontalAlignment=null,
        onBlur,
        onChange,
        onClick,
        onContextMenu,
        onDblClick,
        onMouseDown,
        onMouseEnter,
        onMouseLeave,
        onMouseMove,
        onMouseOut,
        onMouseOver,
        onMouseUp,
        ...otherArgs
    }) {
        let justify, align, text, theClass
        
        if (!row) {
            if (center) {
                verticalAlignment = verticalAlignment || "center"
                horizontalAlignment = horizontalAlignment || "center"
            } else {
                verticalAlignment = verticalAlignment || "top"
                horizontalAlignment = horizontalAlignment || "left"
            }
            justify = verticalAlignment
            align = horizontalAlignment
            text = horizontalAlignment
            theClass = classIds.column
        } else {
            if (center) {
                verticalAlignment = verticalAlignment || "center"
                horizontalAlignment = horizontalAlignment || "center"
            } else {
                verticalAlignment = verticalAlignment || "top"
                horizontalAlignment = horizontalAlignment || "left"
            }
            justify = horizontalAlignment
            align = verticalAlignment
            text = horizontalAlignment
            theClass = classIds.row
        }
        
        const element = <button
            class={combineClasses(theClass, otherArgs.class)}
            style={`justify-content: ${translateAlignment(justify)}; align-items: ${translateAlignment(align)}; text-align: ${text}; ${css(style)}; ${css(otherArgs)};`}
            onBlur={onBlur}
            onChange={onChange}
            onClick={onClick}
            onContextMenu={onContextMenu}
            onDblClick={onDblClick}
            onMouseDown={onMouseDown}
            onMouseEnter={onMouseEnter}
            onMouseLeave={onMouseLeave}
            onMouseMove={onMouseMove}
            onMouseOut={onMouseOut}
            onMouseOver={onMouseOver}
            onMouseUp={onMouseUp}
            {...otherArgs}
            >
                {children}
        </button>

        hoverStyleHelper({ element, hoverStyle })
        
        return element
}