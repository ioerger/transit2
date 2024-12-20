#!/usr/bin/env sh
"\"",`$(echo --% ' |out-null)" >$null;function :{};function dv{<#${/*'>/dev/null )` 2>/dev/null;dv() { #>
echo "1.36.1"; : --% ' |out-null <#'; }; version="$(dv)"; deno="$HOME/.deno/$version/bin/deno"; if [ -x "$deno" ]; then  exec "$deno" run -q -A "$0" "$@";  elif [ -f "$deno" ]; then  chmod +x "$deno" && exec "$deno" run -q -A "$0" "$@";  fi; bin_dir="$HOME/.deno/$version/bin"; exe="$bin_dir/deno"; has () { command -v "$1" >/dev/null; } ;  if ! has unzip; then if ! has apt-get; then  has brew && brew install unzip; else  if [ "$(whoami)" = "root" ]; then  apt-get install unzip -y; elif has sudo; then  echo "Can I install unzip for you? (its required for this command to work) ";read ANSWER;echo;  if [ "$ANSWER" =~ ^[Yy] ]; then  sudo apt-get install unzip -y; fi; elif has doas; then  echo "Can I install unzip for you? (its required for this command to work) ";read ANSWER;echo;  if [ "$ANSWER" =~ ^[Yy] ]; then  doas apt-get install unzip -y; fi; fi;  fi;  fi;  if ! has unzip; then  echo ""; echo "So I couldn't find an 'unzip' command"; echo "And I tried to auto install it, but it seems that failed"; echo "(This script needs unzip and either curl or wget)"; echo "Please install the unzip command manually then re-run this script"; exit 1;  fi;  repo="denoland/deno"; if [ "$OS" = "Windows_NT" ]; then target="x86_64-pc-windows-msvc"; else :;  case $(uname -sm) in "Darwin x86_64") target="x86_64-apple-darwin" ;; "Darwin arm64") target="aarch64-apple-darwin" ;; "Linux aarch64") repo="LukeChannings/deno-arm64" target="linux-arm64" ;; "Linux armhf") echo "deno sadly doesn't support 32-bit ARM. Please check your hardware and possibly install a 64-bit operating system." exit 1 ;; *) target="x86_64-unknown-linux-gnu" ;; esac; fi; deno_uri="https://github.com/$repo/releases/download/v$version/deno-$target.zip"; exe="$bin_dir/deno"; if [ ! -d "$bin_dir" ]; then mkdir -p "$bin_dir"; fi;  if ! curl --fail --location --progress-bar --output "$exe.zip" "$deno_uri"; then if ! wget --output-document="$exe.zip" "$deno_uri"; then echo "Howdy! I looked for the 'curl' and for 'wget' commands but I didn't see either of them. Please install one of them, otherwise I have no way to install the missing deno version needed to run this code"; exit 1; fi; fi; unzip -d "$bin_dir" -o "$exe.zip"; chmod +x "$exe"; rm "$exe.zip"; exec "$deno" run -q -A "$0" "$@"; #>}; $DenoInstall = "${HOME}/.deno/$(dv)"; $BinDir = "$DenoInstall/bin"; $DenoExe = "$BinDir/deno.exe"; if (-not(Test-Path -Path "$DenoExe" -PathType Leaf)) { $DenoZip = "$BinDir/deno.zip"; $DenoUri = "https://github.com/denoland/deno/releases/download/v$(dv)/deno-x86_64-pc-windows-msvc.zip";  [Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12;  if (!(Test-Path $BinDir)) { New-Item $BinDir -ItemType Directory | Out-Null; };  Function Test-CommandExists { Param ($command); $oldPreference = $ErrorActionPreference; $ErrorActionPreference = "stop"; try {if(Get-Command "$command"){RETURN $true}} Catch {Write-Host "$command does not exist"; RETURN $false}; Finally {$ErrorActionPreference=$oldPreference}; };  if (Test-CommandExists curl) { curl -Lo $DenoZip $DenoUri; } else { curl.exe -Lo $DenoZip $DenoUri; };  if (Test-CommandExists curl) { tar xf $DenoZip -C $BinDir; } else { tar -Lo $DenoZip $DenoUri; };  Remove-Item $DenoZip;  $User = [EnvironmentVariableTarget]::User; $Path = [Environment]::GetEnvironmentVariable('Path', $User); if (!(";$Path;".ToLower() -like "*;$BinDir;*".ToLower())) { [Environment]::SetEnvironmentVariable('Path', "$Path;$BinDir", $User); $Env:Path += ";$BinDir"; } }; & "$DenoExe" run -q -A "$PSCommandPath" @args; Exit $LastExitCode; <# 
# */0}`;
import { FileSystem } from "https://deno.land/x/quickr@0.6.56/main/file_system.js"
import { Console, bold, lightRed, yellow, gray } from "https://deno.land/x/quickr@0.6.56/main/console.js"
import { run, Timeout, Env, Cwd, Stdin, Stdout, Stderr, Out, Overwrite, AppendTo, throwIfFails, returnAsString, zipInto, mergeInto } from "https://deno.land/x/quickr@0.6.56/main/run.js"
// import { run, Timeout, Env, Cwd, Stdin, Stdout, Stderr, Out, Overwrite, AppendTo, throwIfFails, returnAsString, zipInto, mergeInto } from "/Users/jeffhykin/repos/quickr/main/run.js"
import { recursivelyOwnKeysOf, get, set, hasKeyList, hasDirectKeyList, remove, merge, compareProperty, recursivelyIterateOwnKeysOf } from "https://deno.land/x/good@1.5.1.0/object.js"

const argsWereGiven = Deno.args.length > 0

const testsFolder = `${FileSystem.thisFolder}/../tests`
const logsFolder = `${FileSystem.thisFolder}/../logs`

import { parse } from "https://deno.land/std@0.168.0/flags/mod.ts"
const flags = parse(Deno.args, {
    boolean: ["dontSave", "cliOnly"],
    string: ["exclude"],
    default: {
        exclude: "",
        cliOnly: false,
    },
})

const runASpecificTest = flags._.length > 0
if (!runASpecificTest && !flags.cliOnly) {
    let oneFailed = false
    for (let pythonTest of ["test_analysis_methods.py","test_norm_methods.py","test_pytransit_tools.py","test_tpp.py"]) {
        console.log('#')
        console.log(`# ${pythonTest}`)
        console.log('#')
        var { success } = await run`python3 ${`./tests/${pythonTest}`}`
        if (!success) {
            oneFailed = true
        }
    }

    if (oneFailed) {
        console.log( "#")
        console.log( "# one of the tests above failed")
        console.log( "#")
    }
}

const paths = (await FileSystem.listItemsIn("tests/cli_tests/")).sort(
    compareProperty({keyList:["basename"]})
)
for (const eachItem of paths) {
    if (runASpecificTest && FileSystem.basename(eachItem.path) != flags._[0]) {
        continue
    }
    if (flags.exclude.split(",").includes(eachItem.basename)) {
        continue
    }
    // clear outputs and results
    for (const eachFile of await FileSystem.listFileItemsIn(eachItem.path)) {
        if (eachFile.basename.match(/\.output|.+\.result($|\.)/)) {
            await FileSystem.remove(eachFile.path)
        }
    }
    // run all the tests
    for (const eachFile of await FileSystem.listFileItemsIn(eachItem.path)) {
        try {
            const permissions = await FileSystem.getPermissions({path:eachFile.path})
            const isExecutable = (permissions.owner.canExecute || permissions.group.canExecute || permissions.others.canExecute)
            if (!isExecutable) {
                continue
            }
            let outputPath
            if (!flags.dontSave) {
                outputPath = `${eachFile.path}.output`
            } else {
                outputPath = `${eachFile.path}.ignore.output`
            }
            await Deno.stdout.write(new TextEncoder().encode(`    Running: ${JSON.stringify(outputPath)}\r`))
            var {success} = await run`${eachFile.path} ${Out(Overwrite(outputPath))}`
            if (success) {
                console.log(`    Passed: ${JSON.stringify(outputPath)}  `)
            } else {
                console.log(`    FAILED: ${JSON.stringify(outputPath)}  `)
            }
        } catch (error) {
            console.error(`    Error running: ${eachFile}`)
        }
    }
}
// (this comment is part of deno-guillotine, dont remove) #>